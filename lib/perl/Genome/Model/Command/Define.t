#!/usr/bin/env genome-perl

use strict;
use warnings;

use above "Genome";
use Test::More;
use Genome::Test::Factory::ProcessingProfile::ReferenceAlignment;

$ENV{UR_DBI_NO_COMMIT} = 1;
$ENV{UR_USE_DUMMY_AUTOGENERATED_IDS} = 1;

my $default_subject_name = 'H_GV-933124G-S.9017';
my $processing_profile = Genome::Test::Factory::ProcessingProfile::ReferenceAlignment->setup_object;

my $tmp_dir = File::Temp::tempdir(CLEANUP => 1);

# test create for a genome model with defined model_name
test_model_from_params(
    model_params => {
        model_name              => "test_model_1",
        subject_name            => $default_subject_name,
        processing_profile      => $processing_profile,
        reference_sequence_build => '93636924', #NCBI-human build 36
    },
);

# test create with a different reference_sequence_build name
test_model_from_params(
    model_params => {
        model_name              => "test_model_2",
        subject_name            => $default_subject_name,
        processing_profile      => $processing_profile,
        reference_sequence_build => '102835775'
    },
);

my $group1 = Genome::ModelGroup->create(name => "test 1");
my $group2 = Genome::ModelGroup->create(name => "test 2");

my $groups = join ",", ($group1->name, $group2->name);
test_model_from_params_with_group($groups);

done_testing();


########################################################3

my $cnt = 0;
sub test_model_from_params_with_group {
    my $model_group_name_string = shift;
    # Get all convergence models associated with each model group we are testing and turn off their automatic building 
    # so nothing happens when we add to a model group in the next test
    my @model_group_names = split ",", $model_group_name_string;
    my @groups;
    for my $model_group_name (@model_group_names) {
        push @groups, Genome::ModelGroup->get(name => $model_group_name); 
        my @convergence_models = Genome::Model::Convergence->get(name => $model_group_name);
        for my $convergence_model (@convergence_models) {
            $convergence_model->auto_build_alignments(0);
        }
    }

    # test normal model and processing profile creation for reference alignment with a model group addition
    test_model_from_params(
        model_params => {
            model_name              => "test_model_incomplete_data_dir_" . Genome::Sys->username,
            subject_name            => $default_subject_name,
            processing_profile      => $processing_profile,
            reference_sequence_build => '93636924', #NCBI-human build 36
            groups => \@groups,
        },
    );
}


sub test_model_from_params {
    my %params = @_;
    my %test_params = %{$params{'test_params'}} if defined $params{'test_params'};

    diag("Test: ".++$cnt);
    my %model_params = %{$params{'model_params'}};
    if ($test_params{'fail'}) {
        &failed_create_model($test_params{'fail'},\%model_params);
    } else {
        &successful_create_model(\%model_params);
    }
}

sub successful_create_model {
    my $params = shift;
    my %params = %{$params};

    my $pp = $params{processing_profile};
    isa_ok($pp,'Genome::ProcessingProfile');

    my $subclass = join('', map { ucfirst($_) } split('\s+',$pp->type_name));
    if (!$params{subject_name}) {
        $params{subject_name} = 'invalid_subject_name';
    }
    my $expected_user_name = Genome::Sys->username;
    my $current_time = UR::Context->current->now;
    my ($expected_date) = split('\w',$current_time);
  
    my $define_class = $pp->class;
    $define_class =~ s/Genome::ProcessingProfile:://;
    $define_class =~ s/::.*//;
    $define_class = "Genome::Model::Command::Define::$define_class";

    my $create_command = $define_class->create(%params);
    isa_ok($create_command,'Genome::Model::Command::Define::Helper');

    $create_command->dump_error_messages(0);
    $create_command->dump_warning_messages(0);
    $create_command->dump_status_messages(0);
    $create_command->queue_error_messages(1);
    $create_command->queue_warning_messages(1);
    $create_command->queue_status_messages(1);

    ok($create_command->execute, 'create command execution successful');
    my @error_messages = $create_command->error_messages();
    print @error_messages, "\n";
    my @warning_messages = $create_command->warning_messages();
    my @status_messages = $create_command->status_messages();
    ok(! scalar(@error_messages), 'no error messages');
    if ($params{'model_name'}) {
        SKIP: {
            skip '_build_model_filesystem paths got moved into Genome::Model', 2 unless (0);
            ok(scalar(@warning_messages), 'create model generated a warning message');
        like($warning_messages[0], qr(model symlink.*already exists), 'Warning message complains about the model link already existing');
        }
    } else {
        ok(!scalar(grep { not m/already exists/ } @warning_messages), 'no warning messages');
        if (@warning_messages) {
            print join("\n",@warning_messages);
        }
    }
    ok(scalar(@status_messages), 'There was a status message');

    my @create_status_messages = grep { /Created model:/ } @status_messages;
    ok(@create_status_messages, 'Got create status message');
    # FIXME - some of those have a second message about creating a directory
    # should probably test for that too
    delete($params{bare_args});
    delete($params{model_name});
    delete($params{reference_sequence_build}); #This property will be the build, not the name/ID
    my $model_id = $create_command->result_model_id;
    ok($model_id, 'got created model id') or die;
    my $model = Genome::Model->get($model_id,);
    isa_ok($model,'Genome::Model::'. $subclass);
    ok($model, 'creation worked for '. $model->name .' model');
    for my $property_name (keys %params) {
        # Don't test this one, since it comes in as a string and gets split. They will not be equal
        next if ($property_name eq "groups");
        is($model->$property_name,$params{$property_name},$property_name .' model indirect accessor');
    }
    is($model->run_as,$expected_user_name,'model run_as accesssor');
    is($model->created_by,$expected_user_name,'model created_by accesssor');
    like($model->creation_date,qr/$expected_date/,'model creation_date accessor');
    is($model->processing_profile_id,$pp->id,'model processing_profile_id indirect accessor');
    is($model->type_name,$pp->type_name,'model type_name indirect accessor');


    # test that model group membership is as expected
    SKIP: {
        skip 'only test group membership if one is expected', 1 unless $params{groups};
        my @groups_expected = @{$params{groups}};
        my @groups_actual = $model->model_groups;
        is(scalar(@groups_actual), scalar(@groups_expected), "Model is a member of the correct number of groups");
    }

    for my $param ($pp->params) {
        my $accessor = $param->name;
        my $value = $pp->$accessor;
        if ($accessor eq 'read_aligner_name' && $value =~ /^maq/) {
            $value = 'maq';
        }
        is($model->$accessor,$value,$accessor .' model indirect accessor');
    }

}

sub failed_create_model {
    my $reason = shift;
    my $params = shift;
    my %params = %{$params};
    my  $create_command = Genome::Model::Command::Define::ReferenceAlignment->create(%params);
    isa_ok($create_command,'Genome::Model::Command::Define::Helper');

    $create_command->dump_error_messages(0);
    $create_command->dump_warning_messages(0);
    $create_command->dump_status_messages(0);
    $create_command->dump_usage_messages(0);
    $create_command->queue_error_messages(1);
    $create_command->queue_warning_messages(1);
    $create_command->queue_status_messages(1);
    {
        *OLD = *STDOUT;
        my $variable;
        open OUT ,'>',\$variable;
        *STDOUT = *OUT;
        ok(!$create_command->execute, 'create command execution failed');
        *STDOUT = *OLD;
    };

    my @error_messages = $create_command->error_messages();
    my @warning_messages = $create_command->warning_messages();
    my @status_messages = $create_command->status_messages();
    ok(scalar(@error_messages), 'There are error messages');
    #like($error_messages[0], qr($reason), 'Error message about '. $reason);
    ok(!scalar(@warning_messages), 'no warning message');
    ok(!scalar(@status_messages), 'no status message');
}

1;

