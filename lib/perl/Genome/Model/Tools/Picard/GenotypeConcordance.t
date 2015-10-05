#!/usr/bin/env genome-perl

BEGIN {
    $ENV{UR_DBI_NO_COMMIT} = 1;
    $ENV{UR_USE_DUMMY_AUTOGENERATED_IDS} = 1;
    $ENV{NO_LSF} = 1;
}

use strict;
use warnings;

use above "Genome";

use Test::More tests => 8;
use Genome::Utility::Test 'compare_ok';

my $pkg = "Genome::Model::Tools::Picard::GenotypeConcordance";
use_ok($pkg);

my $test_dir = __FILE__.".d";

# Params
my $sample_name = 'H_NJ-HCC1395-HCC1395';
my $min_dp = '4';
my $picard_version = '1.123';

# Expected Outputs
my $expected_basename = $test_dir .'/'. $sample_name .'-min_dp'. $min_dp .'-v'. $picard_version;
my $expected_detailed_metrics = $expected_basename .'.genotype_concordance_detail_metrics';
my $expected_summary_metrics = $expected_basename .'.genotype_concordance_summary_metrics';

# Inputs

# The original files were downsampled at this rate
my $size = '0.01';

my $truth_vcf = $test_dir .'/'. $sample_name .'-microarray-'. $size .'.vcf';
my $call_vcf = $test_dir .'/'. $sample_name .'-calls_x_microarray-'. $size .'.vcf';

# Outputs
my $output = Genome::Sys->create_temp_file_path($sample_name);    

# Create test command
my $cmd = Genome::Model::Tools::Picard::GenotypeConcordance->create(
    truth_vcf => $truth_vcf,
    call_vcf => $call_vcf,
    output => $output,
    truth_sample => $sample_name,
    call_sample => $sample_name,
    min_dp => $min_dp,
    use_version => $picard_version,
);
ok($cmd,'create GenotypeConcorndance command');

# Execute test command
ok($cmd->execute,'execute '. $sample_name .' GenotypeConcordance');

# Validate output files were created
my $output_detailed_metrics = $output .'.genotype_concordance_detail_metrics';
my $output_summary_metrics = $output .'.genotype_concordance_summary_metrics';
ok(-e $output_summary_metrics, 'summary metrics exist');
ok(-e $output_detailed_metrics, 'detailed metrics exist');       


# Compare output files to expected
my $filters = [
   qr(^# .*$),
];
compare_ok($expected_detailed_metrics,$output_detailed_metrics, name => 'detailed metrics files match', filters => $filters );
compare_ok($expected_summary_metrics,$output_summary_metrics, name => 'summary metrics files match', filters => $filters );

# test the parsing of the metrics
my $expected_summary_metrics_hashref = Genome::Model::Tools::Picard::GenotypeConcordance->parse_file_into_metrics_hashref($expected_summary_metrics);
my $output_summary_metrics_hashref = Genome::Model::Tools::Picard::GenotypeConcordance->parse_file_into_metrics_hashref($output_summary_metrics);
is_deeply($output_summary_metrics_hashref, $expected_summary_metrics_hashref, 'expected and output summary metrics are the same');
