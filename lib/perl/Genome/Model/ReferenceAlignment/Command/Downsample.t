#!/usr/bin/env genome-perl

use strict;
use warnings;

BEGIN {
    $ENV{UR_DBI_NO_COMMIT} = 1;
    $ENV{UR_USE_DUMMY_AUTOGENERATED_IDS} = 1;
    $ENV{NO_LSF} = 1;
};

use above 'Genome';
use Test::More;

use Genome::Test::Factory::InstrumentData::MergedAlignmentResult;
use Genome::Test::Factory::Model::ReferenceAlignment;
use Genome::Test::Factory::InstrumentData::Solexa;
use Genome::Test::Factory::Library;
use Genome::Test::Factory::Sample;
use Genome::Test::Factory::Build;

my $pkg = 'Genome::Model::ReferenceAlignment::Command::Downsample';
use_ok($pkg);

my $test_dir = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-ReferenceAlignment-Command-Downsample';

my $sample  = Genome::Test::Factory::Sample->setup_object(name => 'test_sample_1', source_id => $pkg);
my $library = Genome::Test::Factory::Library->setup_object(sample => $sample);
my $instrument_data = Genome::Test::Factory::InstrumentData::Solexa->setup_object(library => $library);

my $model = Genome::Test::Factory::Model::ReferenceAlignment->setup_object(
    instrument_data => [$instrument_data],
);

my $build = Genome::Test::Factory::Build->setup_object(
    model_id => $model->id,
    status => 'Succeeded',
);

my $merged_alignment_result = Genome::Test::Factory::InstrumentData::MergedAlignmentResult->setup_object(
    output_dir => $test_dir.'/alignments',
    id => '102922275_merged_rmdup',
);

Sub::Install::reinstall_sub({
    into => "Genome::Model::Build::ReferenceAlignment",
    as   => 'merged_alignment_result',
    code => sub {
        return $merged_alignment_result;
    },
});

my $bam_file = $test_dir.'/alignments/102922275_merged_rmdup.bam';
is($build->whole_rmdup_bam_file, $bam_file,  "bam path correct");
ok($build->isa("Genome::Model::Build::ReferenceAlignment"), "Generated a build");

my $downsample_ratio = 0.5;

my $cmd = Genome::Model::ReferenceAlignment::Command::Downsample->create(
    model => $model,
    coverage_in_ratio => $downsample_ratio,
);
isa_ok($cmd, 'Genome::Model::ReferenceAlignment::Command::Downsample', 'created command');
ok($cmd->execute, 'executed command');

my $new_model = $cmd->_new_model;
my $new_model_name = $model->name . '_downsample_' . $downsample_ratio;
isa_ok($new_model, 'Genome::Model::ReferenceAlignment', 'the command produced a new model as expected');
is($new_model->name, $new_model_name, "New model is named correctly: $new_model_name");

my @instrument_data = $new_model->instrument_data;
my $new_lib_name = $library->name.'-extlibs';
is($instrument_data[0]->library->name, $new_lib_name, "Imported instrument data library is named ok: $new_lib_name");

done_testing();
