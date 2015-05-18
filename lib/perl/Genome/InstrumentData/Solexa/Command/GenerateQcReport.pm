package Genome::InstrumentData::Solexa::Command::GenerateQcReport;

use strict;
use warnings;

use Genome;

class Genome::InstrumentData::Solexa::Command::GenerateQcReport {
    is => 'Command::V2',
    doc => '', #TODO: write me
};

sub execute {
    my $self = shift;
    # XXX The query and next statement below are heuristics to limit this to only X Ten Data.
    my $iter = Genome::InstrumentData::Solexa->create_iterator(read_length => 151);
    print $self->header;
    while(my $inst_data = $iter->next) {
        my $guard = UR::Context::AutoUnloadPool->create();
        next unless $inst_data->run_name =~ /E00/;

        my %attributes_hash = create_attributes_hash($inst_data);
        my ($alignment_result) = $inst_data->get_default_alignment_results; #only grabbing one.
        next unless $alignment_result;
        my $db_metrics_hash = create_metrics_hash($alignment_result) if $alignment_result;
        my ($bamqc_result) = Genome::InstrumentData::AlignmentResult::Merged::BamQc->get(alignment_result_id => $alignment_result->id);
        next unless $bamqc_result;
        my $bamqc_metrics = create_metrics_hash($bamqc_result) if $bamqc_result;
        print join("\t",
            #info about where this was run
            $inst_data->id,
            $self->_illumina_machine_serial_number($inst_data->run_name),
            $inst_data->flow_cell_id,
            $inst_data->lane,
            $inst_data->run_start_date_formatted,

            #what was run
            $inst_data->index_sequence ? $inst_data->index_sequence : 'unknown',
            $inst_data->analysis_projects ? $inst_data->analysis_projects->name : '',
            $inst_data->sample_name,
            $inst_data->library_name,
            # TODO consider splitting out read1 and read2 as separate lines. Might make data mining easier.
            #data stats, base qualities, yield, error rates
            $inst_data->fwd_read_length,
            $inst_data->fwd_clusters,
            $inst_data->fwd_kilobases_read,
            _calculate_avg_qscore(\%attributes_hash, 'fwd', $inst_data->fwd_clusters),
            _calculate_percent_gt_q30(\%attributes_hash, 'fwd', $inst_data->fwd_clusters),
            $inst_data->read_1_pct_mismatch,
            $inst_data->rev_read_length,
            $inst_data->rev_clusters,
            $inst_data->rev_kilobases_read,
            _calculate_avg_qscore(\%attributes_hash, 'rev', $inst_data->rev_clusters),
            _calculate_percent_gt_q30(\%attributes_hash, 'rev', $inst_data->rev_clusters),
            $inst_data->read_2_pct_mismatch,

            #alignments
            $bamqc_metrics->{'bam_qc-FlagstatMetrics-reads_mapped_percentage'},
            $bamqc_metrics->{'bam_qc-FlagstatMetrics-reads_mapped_in_proper_pairs_percentage'},
            $bamqc_metrics->{'bam_qc-FlagstatMetrics-reads_mapped_as_singleton_percentage'},
            $bamqc_metrics->{'bam_qc-FlagstatMetrics-reads_mapped_in_pair'} ? sprintf("%0.02f",$bamqc_metrics->{'bam_qc-FlagstatMetrics-hq_reads_mapped_in_interchromosomal_pairs'} / $bamqc_metrics->{'bam_qc-FlagstatMetrics-reads_mapped_in_pair'} * 100) : 0,
            sprintf("%0.02f",$bamqc_metrics->{'bam_qc-AlignmentSummaryMetrics-CATEGORY-PAIR-PCT_PF_READS_ALIGNED'} * 100),
            sprintf("%0.02f",$bamqc_metrics->{'bam_qc-AlignmentSummaryMetrics-CATEGORY-PAIR-PCT_CHIMERAS'} * 100),
            sprintf("%0.02f",$bamqc_metrics->{'bam_qc-AlignmentSummaryMetrics-CATEGORY-PAIR-PF_HQ_ERROR_RATE'} * 100),
            sprintf("%0.02f",$bamqc_metrics->{'bam_qc-AlignmentSummaryMetrics-CATEGORY-PAIR-PF_INDEL_RATE'} * 100),
            sprintf("%0.02f",$bamqc_metrics->{'bam_qc-AlignmentSummaryMetrics-CATEGORY-PAIR-PF_MISMATCH_RATE'} * 100),
            sprintf("%0.02f",$bamqc_metrics->{'bam_qc-AlignmentSummaryMetrics-CATEGORY-PAIR-STRAND_BALANCE'} * 100),


            #library insert size
            $bamqc_metrics->{'bam_qc-InsertSizeMetrics-PAIR_ORIENTATION-FR-MEDIAN_INSERT_SIZE'},
            $bamqc_metrics->{'bam_qc-InsertSizeMetrics-PAIR_ORIENTATION-FR-MEDIAN_ABSOLUTE_DEVIATION'},
            $bamqc_metrics->{'bam_qc-InsertSizeMetrics-PAIR_ORIENTATION-FR-MEAN_INSERT_SIZE'},
            $bamqc_metrics->{'bam_qc-InsertSizeMetrics-PAIR_ORIENTATION-FR-STANDARD_DEVIATION'},

            #gc bias
            $bamqc_metrics->{'bam_qc-GcBiasSummary-WINDOW_SIZE-100-AT_DROPOUT'},
            $bamqc_metrics->{'bam_qc-GcBiasSummary-WINDOW_SIZE-100-GC_DROPOUT'},

        ), "\n";
        UR::Context->clear_cache(dont_unload => ['Genome::InstrumentData::Solexa']);
    }
    return 1;
}

sub _illumina_machine_serial_number {
    my ($self, $run_name) = @_;
    # XXX This could be a method on InstrumentData::Solexa if we want.
    my ($date, $machine, $run_id, $flow_cell) = split "_", $run_name;
    unless(defined $date
        && defined $machine
        && defined $run_id
        && defined $flow_cell) {
        die "Unable to parse run_name properly from $run_name\n";
    }
    return $machine;
}

sub create_metrics_hash {
    my $object_with_metrics = shift;
    # TODO Check if we can even call metrics.
    my %metrics;
    map { $metrics{$_->metric_name} = $_->metric_value } $object_with_metrics->metrics;
    return \%metrics;
}

sub create_attributes_hash {
    my $instrument_data = shift;
    my @attributes = $instrument_data->attributes;;
    my %attributes_hash;
    for my $a (@attributes){
        $attributes_hash{$a->attribute_label} = $a->attribute_value;
    }
    return %attributes_hash;
}

sub _calculate_avg_qscore {
    my ($attributes_hash, $prefix, $clusters) = @_;
    my $base_quality_sum = $attributes_hash->{$prefix . '_base_quality_sum'};
    my $read_length = $attributes_hash->{'fwd_read_length'};
    if($base_quality_sum && $read_length) {
        return (($base_quality_sum / ($clusters * $read_length)) / 100);
    }
    else {
        return '';
    }
}

sub _calculate_percent_gt_q30 {
    my ($attributes_hash, $prefix, $clusters) = @_;
    my $q30_base_count = $attributes_hash->{$prefix . '_q30_base_count'};
    my $read_length = $attributes_hash->{'fwd_read_length'};
    if($q30_base_count && $read_length) {
        return ($q30_base_count / ($clusters * $read_length));
    }
    else {
        return '';
    }
}

sub header {
        return join("\t",
            #info about where this was run
            "instrument_data_id",
            "machine",
            "flow_cell_id",
            "lane",
            "run_date",
            "index_sequence",
            "analysis_project",
            "sample_name",
            "library_name",
            "fwd_read_length",
            "fwd_clusters",
            "fwd_kilobases_read",
            "fwd_avg_base_qscore",
            "fwd_percent_bases_gt_q30",
            "fwd_mismatch_rate",
            "rev_read_length",
            "rev_clusters",
            "rev_kilobases_read",
            "rev_avg_base_qscore",
            "rev_percent_bases_gt_q30",
            "rev_mismatch_rate",
            "flagstat_map_pct",
            "flagstat_map_as_proper_pair_pct",
            "flagstat_map_as_singleton_pct",
            "flagstat_map_interchromosomal_pct",
            #alignments
            "alignment_summary_map_pct",
            "alignment_summary_pct_chimeras",
            "alignment_summary_hq_indel_rate",
            "alignment_summary_mismatch_rate",
            "alignment_summary_strand_balance",
            "median_insert_size",
            "median_absolute_deviation_insert_size",
            "mean_insert_size",
            "standard_deviation_insert_size",
            "at_dropout",
            "gc_dropout",
        ). "\n";
}
    
1;
