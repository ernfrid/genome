#!/usr/bin/env genome-perl

BEGIN {
    $ENV{UR_DBI_NO_COMMIT} = 1;
    $ENV{UR_USE_DUMMY_AUTOGENERATED_IDS} = 1;
}

use strict;
use warnings;

use above "Genome";
use Test::More;
use Sub::Override;
use Genome::Test::Factory::InstrumentData::Solexa;
use Genome::Test::Factory::InstrumentData::AlignmentResult;
use Cwd qw(abs_path);

my $pkg = 'Genome::Qc::Tool::VerifyBamId';
use_ok($pkg);

my $data_dir = __FILE__.".d";

use Genome::Qc::Tool;
my $sample_name_override = Sub::Override->new(
    'Genome::Qc::Tool::sample_name',
    sub { return 'TEST-patient1-somval_tumor1'; },
);

my $instrument_data = Genome::Test::Factory::InstrumentData::Solexa->setup_object(
    flow_cell_id => '12345ABXX',
    lane => '2',
    subset_name => '2',
    run_name => 'example',
    id => 'NA12878',
);
my $alignment_result = Genome::Test::Factory::InstrumentData::AlignmentResult->setup_object(
    instrument_data => $instrument_data,
);

my $vcf_file = abs_path(File::Spec->join($data_dir, 'Omni25_genotypes_1525_samples_v2.b37.PASS.ALL.sites.chrY.vcf'));
my $bam_file = abs_path(File::Spec->join($data_dir, 'speedseq_merged.bam'));

use Genome::Qc::Config;
my $config_override = Sub::Override->new(
    'Genome::Qc::Config::get_commands_for_alignment_result',
    sub {
        return {
            verify_bam_id => {
                class => 'Genome::Qc::Tool::VerifyBamId',
                params => {
                    vcf => $vcf_file,
                    bam => $bam_file,
                    max_depth => '150',
                    precise => '1',
                    version => '20120620',
                    ignore_read_group => 0,
                },
            }
        };
    },
);

my $command = Genome::Qc::Run->create(
    config_name => 'testing-qc-run',
    alignment_result => $alignment_result,
    %{Genome::Test::Factory::SoftwareResult::User->setup_user_hash},
);
ok($command->execute, "Command executes ok");

my %tools = $command->output_result->_tools;
my ($tool) = values %tools;
ok($tool->isa($pkg), 'Tool created successfully');

my $output = $tool->qc_metrics_file;
my @expected_cmd_line = (
    '/usr/bin/verifyBamID20120620',
    '--vcf',
    $vcf_file,
    '--bam',
    $bam_file,
    '--out',
    $output,
    '--maxDepth',
    150,
    '--precise',
);
is_deeply([$tool->cmd_line], [@expected_cmd_line], 'Command line list as expected');

my %expected_metrics = (
    '2883581792-2883255521-#READS' => 0,
    '2883581792-2883255521-#SNPS' => 1172,
    '2883581792-2883255521-AVG_DP' => '0.00',
    '2883581792-2883255521-CHIPLK0' => 'NA',
    '2883581792-2883255521-CHIPLK1' => 'NA',
    '2883581792-2883255521-CHIPMIX' => 'NA',
    '2883581792-2883255521-CHIP_ID' => 'NA',
    '2883581792-2883255521-CHIP_RA' => 'NA',
    '2883581792-2883255521-CHIP_RH' => 'NA',
    '2883581792-2883255521-DPREF' => 'NA',
    '2883581792-2883255521-FREELK0' => '0.00',
    '2883581792-2883255521-FREELK1' => '0.00',
    '2883581792-2883255521-FREEMIX' => '0.00000',
    '2883581792-2883255521-FREE_RA' => 'NA',
    '2883581792-2883255521-FREE_RH' => 'NA',
    '2883581792-2883255521-RDPALT' => 'NA',
    '2883581792-2883255521-RDPHET' => 'NA',
    'ALL-#READS' => 0,
    'ALL-#SNPS' => 1172,
    'ALL-AVG_DP' => '0.00',
    'ALL-CHIPLK0' => 'NA',
    'ALL-CHIPLK1' => 'NA',
    'ALL-CHIPMIX' => 'NA',
    'ALL-CHIP_ID' => 'NA',
    'ALL-CHIP_RA' => 'NA',
    'ALL-CHIP_RH' => 'NA',
    'ALL-DPREF' => 'NA',
    'ALL-FREELK0' => '0.00',
    'ALL-FREELK1' => '0.00',
    'ALL-FREEMIX' => '0.00000',
    'ALL-FREE_RA' => 'NA',
    'ALL-FREE_RH' => 'NA',
    'ALL-RDPALT' => 'NA',
    'ALL-RDPHET' => 'NA',
);
is_deeply({$command->output_result->get_metrics}, {%expected_metrics}, 'Parsed metrics as expected');

done_testing;
