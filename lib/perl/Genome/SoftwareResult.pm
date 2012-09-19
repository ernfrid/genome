package Genome::SoftwareResult;

use strict;
use warnings;

use Genome;
use Digest::MD5 qw(md5_hex);
use Cwd;
use File::Basename qw(fileparse);
use Data::Dumper;

class Genome::SoftwareResult {
    is_abstract => 1,
    table_name => 'SOFTWARE_RESULT',
    subclass_description_preprocessor => 'Genome::SoftwareResult::_expand_param_and_input_properties',
    subclassify_by => 'subclass_name',
    id_by => [
        id => { is => 'Number', len => 20 },
    ],
    attributes_have => [
        is_param => { is => 'Boolean', is_optional=>'1' },
        is_input => { is => 'Boolean', is_optional=>'1' },
        is_metric => { is => 'Boolean', is_optional=>'1' }
    ],
    has => [
        module_version      => { is => 'Text', len => 64, column_name => 'VERSION', is_optional => 1 },
        subclass_name       => { is => 'Text', len => 255, column_name => 'CLASS_NAME' },
        inputs_bx           => { is => 'UR::BoolExpr', id_by => 'inputs_id', is_optional => 1 },
        inputs_id           => { is => 'Text', len => 4000, column_name => 'INPUTS_ID', implied_by => 'inputs_bx', is_optional => 1 },
        params_bx           => { is => 'UR::BoolExpr', id_by => 'params_id', is_optional => 1 },
        params_id           => { is => 'Text', len => 4000, column_name => 'PARAMS_ID', implied_by => 'params_bx', is_optional => 1 },
        output_dir          => { is => 'Text', len => 1000, column_name => 'OUTPUTS_PATH', is_optional => 1 },
        test_name           => { is_param => 1, is_delegated => 1, is_mutable => 1, via => 'params', to => 'value_id', where => ['name' => 'test_name'], is => 'Text', doc => 'Assigns a testing tag to the result.  These will not be used in default processing', is_optional => 1 },
        _lock_name          => { is_param => 1, is_optional => 1, is_transient => 1 },
    ],
    has_many_optional => [
        params              => { is => 'Genome::SoftwareResult::Param', reverse_as => 'software_result'},
        inputs              => { is => 'Genome::SoftwareResult::Input', reverse_as => 'software_result'},
        metrics             => { is => 'Genome::SoftwareResult::Metric', reverse_as => 'software_result'},
        users               => { is => 'Genome::SoftwareResult::User', reverse_as => 'software_result'},
        disk_allocations    => { is => 'Genome::Disk::Allocation', reverse_as => 'owner'},
        build_ids           => { via => 'users', to => 'user_id', } # where => ['user_class_name isa' => 'Genome::Model::Build'] },
    ],
    schema_name => 'GMSchema',
    data_source => 'Genome::DataSource::GMSchema',
    doc => 'base class for managed data sets, with database tracking for params, inputs, metrics, and disk',
};

our %LOCKS;

# You must specify enough parameters to uniquely identify an object to get a result.
# If two users specify different sets of parameters that uniquely identify the result,
# they will create different locks.
sub get_with_lock {
    my $class = shift;

    my $params_processed = $class->_gather_params_for_get_or_create($class->_preprocess_params_for_get_or_create(@_));

    my %is_input = %{$params_processed->{inputs}};
    my %is_param = %{$params_processed->{params}};

    # Only try with lock if object does not exist since locking causes
    # a performance hit. It is assumed that if an object is found it is
    # complete. If this is a bad assumption then we need to add a
    # status to SoftwareResults.
    my $lock;
    my @objects = $class->get(%is_input, %is_param);
    unless (@objects) {
        my $subclass = $params_processed->{subclass};
        unless ($lock = $subclass->_lock(%is_input, %is_param)) {
            die "Failed to get a lock for " . Dumper(\%is_input,\%is_param);
        }

        eval {
            @objects = $class->get(%is_input, %is_param);
        };
        my $error = $@;

        if ($error) {
            $class->error_message('Failed in get! ' . $error);
            $class->_release_lock_or_die($lock, "Failed to unlock during get_with_lock.");
            die $class->error_message;
        }
    }

    if (@objects > 1) {
        $class->error_message("Multiple results returned for SoftwareResult::get_with_lock.  To avoid this, call get_with_lock with enough parameters to uniquely identify a SoftwareResult.");
        $class->error_message("Parameters used for the get: " . Data::Dumper::Dumper %is_input . Data::Dumper::Dumper %is_param);
        $class->error_message("Objects gotten: " . Data::Dumper::Dumper @objects);
        $class->_release_lock_or_die($lock, "Failed to unlock during get_with_lock with multiple results.") if $lock;
        die $class->error_message;
    }

    my $result = $objects[0];
    if ($result && $lock) {
        $result->_lock_name($lock);

        $result->status_message("Cleaning up lock $lock...");
        unless ($result->_unlock) {
            $result->error_message("Failed to unlock after getting software result");
            die "Failed to unlock after getting software result";
        }
        $result->status_message("Cleanup completed for lock $lock.");
    } elsif ($lock) {
        $class->_release_lock_or_die($lock, "Failed to unlock after not finding software result.");
    }

    return $result;
}

sub get_or_create {
    my $class = shift;

    my $params_processed = $class->_gather_params_for_get_or_create($class->_preprocess_params_for_get_or_create(@_));
    my %is_input = %{$params_processed->{inputs}};
    my %is_param = %{$params_processed->{params}};
    
    my @objects = $class->get(%is_input, %is_param);

    unless (@objects) {
        @objects = $class->create(@_);
        unless (@objects) {
            # see if the reason we failed was b/c the objects were created while we were locking...
            @objects = $class->get(%is_input, %is_param);
            unless (@objects) {
                $class->error_message("Could not create a $class for params " . Data::Dumper::Dumper(\@_) . " even after trying!");
                die $class->error_message();
            }
        }
    }
   
    if (@objects > 1) {
        return @objects if wantarray;
        my @ids = map { $_->id } @objects;
        die "Multiple matches for $class but get or create was called in scalar context!  Found ids: @ids";
    } else {
        return $objects[0];
    }
}

sub create {
    my $class = shift;

    if ($class eq __PACKAGE__ || $class->__meta__->is_abstract) {
        # this class is abstract, and the super-class re-calls the constructor from the correct subclass
        return $class->SUPER::create(@_);
    }

    my $params_processed = $class->_gather_params_for_get_or_create($class->_preprocess_params_for_get_or_create(@_));
    my %is_input = %{$params_processed->{inputs}};
    my %is_param = %{$params_processed->{params}};

    my @previously_existing = $class->get(%is_input, %is_param);
    if (@previously_existing > 0) {
        $class->error_message("Attempt to create an $class but it looks like we already have one with those params " . Dumper(\@_));
        return;
    }

    my $lock;
    unless ($lock = $class->_lock(%is_input, %is_param)) {
        die "Failed to get a lock for " . Dumper(\%is_input,\%is_param);
    }

    # TODO; if an exception occurs before this is assigned to the object, we'll have a stray lock
    # We need to ensure that we get cleanup on die.

    # we might have had to wait on the lock, in which case someone else was probably creating that entity
    # do a "reload" here to force another trip back to the database to see if a software result was created
    # while we were waiting on the lock.
    (@previously_existing) = UR::Context->current->reload($class,%is_input,%is_param);

    if (@previously_existing > 0) {
        $class->error_message("Attempt to create an $class but it looks like we already have one with those params " . Dumper(\@_));
        $class->_release_lock_or_die($lock, "Failed to release lock in create before committing SoftwareResult.");
        return; 
    }
    
    # We need to update the indirect mutable accessor logic for non-nullable
    # hang-offs to delete the entry instead of setting it to null.  Otherwise
    # we get SOFTWARE_RESULT_PARAM entries with a NULL, and unsavable PARAM_VALUE.
    # also remove empty strings because that's equivalent to a NULL to the database

    # Do the same for inputs (e.g. alignment results have nullable segment values for instrument data, which are treated as inputs)
    my @param_remove = grep { not (defined $is_param{$_}) || $is_param{$_} eq "" } keys %is_param;
    my @input_remove = grep { not (defined $is_input{$_}) || $is_input{$_} eq "" } keys %is_input;
    my $bx = $class->define_boolexpr($class->_preprocess_params_for_get_or_create(@_));
    for my $i (@param_remove, @input_remove) {
        $bx = $bx->remove_filter($i);
    }

    my $self = $class->SUPER::create($bx);
    unless ($self) {
        $class->_release_lock_or_die($lock,"Failed to unlock during create after committing SoftwareResult.");
        return;
    }

    $self->_lock_name($lock);

    my $unlock_callback = sub {
        $self->_unlock;
    };
    $self->create_subscription(method=>'commit', callback=>$unlock_callback);
    $self->create_subscription(method=>'delete', callback=>$unlock_callback);

    if (my $output_dir = $self->output_dir) {
        if (-d $output_dir) {
            my @files = glob("$output_dir/*");
            if (@files) {
                $self->delete;
                die "Found files in output directory $output_dir!:\n\t" 
                    . join("\n\t", @files);
            }
            else {
                $self->status_message("No files in $output_dir.");
            }
        }
        else {
            $self->status_message("Creating output directory $output_dir...");
            eval {
                Genome::Sys->create_directory($output_dir)
            };
            if ($@) {
                $self->delete;
                die $@;
            }
        }
    }

    $self->module_version($self->resolve_module_version) unless defined $self->module_version;
    $self->subclass_name($class);
    return $self;
}

sub _gather_params_for_get_or_create {
    my $class = shift;
    my $bx = UR::BoolExpr->resolve_normalized_rule_for_class_and_params($class, @_);

    my %params = $bx->params_list;
    my %is_input;
    my %is_param;
    my $class_object = $class->__meta__;
    for my $key ($class->property_names) {
        my $meta = $class_object->property_meta_for_name($key);
        if ($meta->{is_input} && exists $params{$key}) {
            $is_input{$key} = $params{$key};
        } elsif ($meta->{is_param} && exists $params{$key}) {
            $is_param{$key} = $params{$key};
        }

    }

    my $inputs_bx = UR::BoolExpr->resolve_normalized_rule_for_class_and_params($class, %is_input);
    my $params_bx = UR::BoolExpr->resolve_normalized_rule_for_class_and_params($class, %is_param);

    my %software_result_params = (#software_version=>$params_bx->value_for('aligner_version'),
        params_id=>$params_bx->id,
        inputs_id=>$inputs_bx->id,
        subclass_name=>$class
    );

    return {
        software_result_params => \%software_result_params,
        subclass => $class,
        inputs=>\%is_input,
        params=>\%is_param,
    };
}

sub _preprocess_params_for_get_or_create {
    my $class = shift;
    if(scalar @_ eq 1) {
        return @_; #don't process a UR::BoolExpr or plain ID
    }

    my %params = @_;

    my $class_object = $class->__meta__;
    for my $key ($class->property_names) {
        my $meta = $class_object->property_meta_for_name($key);

        for my $t ('input', 'param') {
            if ($meta->{'is_' . $t} && $meta->is_many) {
                my $value_list = delete $params{$key};
                if((defined $value_list) && (scalar @$value_list)) {
                    my @values = sort map { Scalar::Util::blessed($_)? $_->id : $_ } @$value_list;
                    my $t_params = $params{$t . 's'} || [];
                    for my $i (0..$#values) {
                        my $value = $values[$i];
                        push @$t_params, {'name' => $key . '-' . $i, 'value_id' => $value};
                    }
                    $params{$t . 's'} = $t_params;

                    $params{$key . '_count'} = scalar @values;
                    $params{$key . '_md5'} = md5_hex( join(':', @values));
                } else {
                    $params{$key . '_count'} = 0;
                    $params{$key . '_md5'} = undef;
                }
            }
        }
    }
    return %params;
}

sub resolve_module_version {
    my $class = shift;
    my $revision = Genome::Sys->snapshot_revision;
    # the revision may be a series of long paths
    # truncate the revision if it is too long to store, but put a * at the end so we can tell this is the case
    my $pmeta = $class->__meta__->property("module_version");
    my $len = $pmeta->data_length - 1;
    if (length($revision) > $len) {
        $revision = substr($revision,0,$len) . '*';
    }
    return $revision;
}

sub _expand_param_and_input_properties {
    my ($class, $desc) = @_;
    for my $t ('input','param','metric') {
        while (my ($prop_name, $prop_desc) = each(%{ $desc->{has} })) {
            if (exists $prop_desc->{'is_'.$t} and $prop_desc->{'is_'.$t}) {
                my $is_many = ($t ne 'metric' and exists $prop_desc->{'is_many'} and $prop_desc->{'is_many'});

                my $name_name;
                if ($t eq 'metric') {
                    $prop_desc->{'to'} = 'metric_value';
                    $name_name = 'metric_name';
                }
                else {
                    # TODO This logic was borrowed in the Model.pm's _resolve_to_for_prop_desc so 
                    # when this is refactored, that should also be updated.
                    if (exists $prop_desc->{'data_type'} and $prop_desc->{'data_type'}) {
                        my $prop_class = UR::Object::Property->_convert_data_type_for_source_class_to_final_class(
                            $prop_desc->{'data_type'},
                            $class
                        );
                        if ($prop_class->isa("UR::Value")) {
                            $prop_desc->{'to'} = 'value_id';
                        } else {
                            $prop_desc->{'to'} = 'value_obj';
                        }
                    } 
                    else {
                        $prop_desc->{'to'} = 'value_id';
                    }
                    $name_name = 'name';    
                }

                $prop_desc->{'is_delegated'} = 1;

                if($is_many) {
                    $prop_desc->{'where'} = [
                        $name_name . ' like' => $prop_name . '-%',
                    ];
                } 
                else {
                    $prop_desc->{'where'} = [
                        $name_name => $prop_name
                    ];
                }
                

                $prop_desc->{'is_mutable'} = 1;
                $prop_desc->{'via'} = $t.'s';

                if($is_many) {
                    my $md5_name = $prop_name . '_md5';
                    unless(exists $desc->{has}{$md5_name}) {
                        my $md5_prop = {};
                        $md5_prop->{'is'} = 'Text';
                        $md5_prop->{'is_param'} = 1;
                        $md5_prop->{'is_delegated'} = 1;
                        $md5_prop->{'via'} = 'params';
                        $md5_prop->{'to'} = 'value_id';
                        $md5_prop->{'where'} = [ 'name' => $md5_name ];
                        $md5_prop->{'doc'} = 'MD5 sum of the sorted list of values for ' . $prop_name;
                        $md5_prop->{'is_mutable'} = 1;
                        $md5_prop->{'is_optional'} = 1;

                        $md5_prop->{'property_name'} = $md5_name;
                        $md5_prop->{'class_name'} = $desc->{class_name};
                        $desc->{has}{$md5_name} = $md5_prop;
                    }

                    my $count_name = $prop_name . '_count';
                    unless(exists $desc->{has}{$count_name}) {
                        my $count_prop = {};
                        $count_prop->{'is'} = 'Number';
                        $count_prop->{'is_param'} = 1;
                        $count_prop->{'is_delegated'} = 1;
                        $count_prop->{'via'} = 'params';
                        $count_prop->{'to'} = 'value_id';
                        $count_prop->{'where'} = [ 'name' => $count_name ];
                        $count_prop->{'doc'} = 'number of values for ' . $prop_name;
                        $count_prop->{'is_mutable'} = 1;
                        $count_prop->{'is_optional'} = 1;

                        $count_prop->{'property_name'} = $count_name;
                        $count_prop->{'class_name'} = $desc->{class_name};
                        $desc->{has}{$count_name} = $count_prop;
                    }
                }
            }
        }
    }
    return $desc;
}

sub delete {
    my $self = shift;

    my $class_name = $self->class;
    my @users = $self->users;
    if (@users) {
        my $name = $self->__display_name__;
        die "Refusing to delete $class_name $name as it still has users:\n\t"
            .join("\n\t", map { $_->user_class . "\t" . $_->user_id } @users);
    }

    my @to_nuke = ($self->params, $self->inputs, $self->metrics); 

    #If we use any other results, unregister ourselves as users
    push @to_nuke, Genome::SoftwareResult::User->get(user_class_name => $class_name, user_id => $self->id);

    for (@to_nuke) {
        unless($_->delete) {
            die "Failed to delete: " . Data::Dumper::Dumper($_);
        }
    }

    #creating an anonymous sub to delete allocations when commit happens
    my $id = $self->id;
    my $upon_delete_callback = sub { 
        print "Now Deleting Allocation with owner_id = $id\n";
        my $allocation = Genome::Disk::Allocation->get(owner_id=>$id, owner_class_name=>$class_name);
        if ($allocation) {
            $allocation->deallocate; 
        }
    };

    #hook our anonymous sub into the commit callback
    $class_name->ghost_class->add_observer(aspect=>'commit', callback=>$upon_delete_callback);
    
    return $self->SUPER::delete(@_); 
}

sub _lock {
    my $class = shift;
    
    my $resource_lock_name = $class->_resolve_lock_name(@_);

    # if we're already locked, just increment the lock count
    $LOCKS{$resource_lock_name} += 1;
    return $resource_lock_name if ($LOCKS{$resource_lock_name} > 1);
   
    my $lock = Genome::Sys->lock_resource(resource_lock => $resource_lock_name, max_try => 2);
    unless ($lock) {
        $class->status_message("This data set is still being processed by its creator.  Waiting for existing data lock...");
        $lock = Genome::Sys->lock_resource(resource_lock => $resource_lock_name, wait_announce_interval => 600);
        unless ($lock) {
            $class->error_message("Failed to get existing data lock!");
            die($class->error_message);
        }
    }

    return $lock;
}

sub _unlock {
    my $self = shift;

    my $resource_lock_name = $self->_lock_name;
    $self->status_message("Cleaning up lock $resource_lock_name...");

    if (!exists $LOCKS{$resource_lock_name})  {
        $self->error_message("Attempt to unlock $resource_lock_name but this was never locked!");
        die $self->error_message;
    }
    $LOCKS{$resource_lock_name} -= 1;

    return 1 if ($LOCKS{$resource_lock_name} >= 1);

    unless (Genome::Sys->unlock_resource(resource_lock=>$resource_lock_name)) {
        $self->error_message("Couldn't unlock $resource_lock_name.  error message was " . $self->error_message);
        die $self->error_message;
    }

    delete $LOCKS{$resource_lock_name};
    $self->status_message("Cleanup completed for lock $resource_lock_name.");
    return 1;
}

sub _resolve_lock_name {
    my $class = shift;
    my $class_string = $class;
    $class_string =~ s/\:/\-/g;

    my $be = UR::BoolExpr->resolve_normalized($class, @_);
    no warnings;
    my $params_and_inputs_list=join "___", $be->params_list;
    # sub out dangerous directory separators
    $params_and_inputs_list =~ s/\//\./g;
    use warnings;
    my $params_and_inputs_list_hash = md5_hex($params_and_inputs_list);

    my $resource_lock_name = $ENV{GENOME_LOCK_DIR} . "/genome/$class_string/" .  $params_and_inputs_list_hash;
}

# override _resolve_lock_name (for testing) to append username and time
# This override is used to prevent lock collisions when tests are being run concurrently on the same machine.
if ($ENV{UR_DBI_NO_COMMIT}) {
    warn 'Overriding Genome::SoftwareResult::_resolve_lock_name since UR_DBI_NO_COMMIT is on.' . "\n";
    my $suffix = Genome::Sys->username . '_' . time;
    my $original_resolve_lock_name_sub = \&Genome::SoftwareResult::_resolve_lock_name;
    *Genome::SoftwareResult::_resolve_lock_name = sub {
        my $lock_name = &$original_resolve_lock_name_sub(@_);
        $lock_name .= "_$suffix" unless $lock_name =~ /$suffix/;
        return $lock_name;
    };
}

sub metric_names {
    my $class = shift;
    my $meta = $class->__meta__;
    my @properties = grep { $_->{is_metric} } $meta->_legacy_properties();
    my @names = map { $_->property_name } @properties;
    return @names;
}

sub metrics_hash {
    my $self = shift;
    my @names = $self->metric_names;
    my %hash = map { $self->name } @names;
    return %hash;
}

sub generate_expected_metrics {
    my $self = shift;
    my @names = @_;
    unless (@names) {
        @names = $self->metric_names;
    }
    
    # pre-load all metrics
    my @existing_metrics = $self->metrics;
    
    for my $name (@names) {
        my $metric = $self->metric(name => $name);
        if ($metric) {
            $self->status_message(
                $self->display_name . " has metric "
                . $metric->name 
                . " with value "
                . $metric->value
            );
            next;
        }
        my $method = "_calculate_$name";
        unless ($self->can($method)) {
            $self->error_message("No method $method found!");
            die $self->error_message;
        }
        $self->status_message(
            $self->display_name . " is generating a value for metric "
            . $metric->name 
            . "..."
        );
        my $value = $self->$method();
        unless (defined($value)) {
            $self->error_message(
                $self->display_name . " has metric "
                . $metric->name 
                . " FAILED TO CALCULATE A DEFINED VALUE"
            );
            next;
        }
        $self->$metric($value);
        $self->status_message(
            $self->display_name . " has metric "
            . $metric->name 
            . " with value "
            . $metric->value
        );
    }
}

sub _available_cpu_count {
    my $self = shift; 

    # Not running on LSF, allow only one CPU
    if (!exists $ENV{LSB_MCPU_HOSTS}) {
        return 1;
    }

    my $mval = $ENV{LSB_MCPU_HOSTS};
    my @c = split /\s+/, $mval;

    if (scalar @c != 2) {
        $self->error_message("LSB_MCPU_HOSTS environment variable doesn't specify just one host and one CPU count. (value is '$mval').  Is the span[hosts=1] value set in your resource request?");
        die $self->error_message;
    }

    if ($mval =~ m/(\.*?) (\d+)/) {
        return $2; 
    } else {
        $self->error_message("Couldn't parse the LSB_MCPU_HOSTS environment variable (value is '$mval'). "); 
        die $self->error_message;
    }
    
}

sub _resolve_param_value_from_text_by_name_or_id {
    my $class = shift;
    my $param_arg = shift;

    #First try default behaviour of looking up by name or id
    my @results = Command::V2->_resolve_param_value_from_text_by_name_or_id($class, $param_arg);

    #If that didn't work, and the argument is a filename, see if it's part of our output directory.
    if(!@results and -f $param_arg) {
        my $abs_path = Cwd::abs_path($param_arg);
        my (undef, $dir) = fileparse($abs_path);
        $dir =~ s!/$!!; #remove trailing slash!
        @results = Genome::SoftwareResult->get(output_dir => $dir);
    }

    return @results;
}

sub _release_lock_or_die {
    my ($class, $lock, $error_message) = @_;

    $class->status_message("Cleaning up lock $lock...");

    unless (Genome::Sys->unlock_resource(resource_lock=>$lock)) {
        $class->error_message($error_message);
        die $error_message;
    }
    delete $LOCKS{$lock};

    $class->status_message("Cleanup completed for lock $lock.");
}

1;
