package Genome::Model::Build::ReferenceSequence::Buckets;

use strict;
use warnings;

use Genome;
use Algorithm::Bucketizer;
use List::Util qw(max);

class Genome::Model::Build::ReferenceSequence::Buckets {
    is => 'Genome::SoftwareResult::StageableSimple',
    has_input => [
        reference_sequence_build => {
            is => 'Genome::Model::Build::ReferenceSequence',
            doc => 'The reference to bucketize',
        },
    ],
    has_metric => [
        count => {
            is => 'Number',
            doc => 'The number of buckets for this reference',
        },
    ],
};

sub bucket_list {
    my $self = shift;

    my @buckets;
    for my $bucket (1..$self->count) {
        push @buckets, $self->bucket($bucket);
    }

    return \@buckets;
}

sub bucket {
    my $self = shift;
    my $bucket = shift;

    my $output_dir = $self->output_dir;

    my $bucket_file = File::Spec->join($output_dir, $bucket);
    my @items = Genome::Sys->read_file($bucket_file);
    chomp @items;

    return \@items;
}

sub _is_grch38_primary {
    my $chr_name = shift;
    # see ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/README.20150309.GRCh38_full_analysis_set_plus_decoy_hla
    if ($chr_name =~ /(_alt|_patch|_decoy|EBV)$|^HLA-/) {
        return 0;
    }
    else {
        return 1;
    }
}

sub _run {
    my $self = shift;

    # The below is a horrible hack to generate buckets for my reference sequence
    my $all_chr_lengths = $self->reference_sequence_build->chromosomes_with_lengths;
    my @chromosomes_we_care_about = grep { _is_grch38_primary($_->[0]) } @$all_chr_lengths;
    my $chr_lengths = \@chromosomes_we_care_about;

    my $max_length = max( map $_->[1], @$chr_lengths );

    my $bucketizer = Algorithm::Bucketizer->new(bucketsize => $max_length);
    for my $chr_with_length (@$chr_lengths) {
        $bucketizer->add_item( @$chr_with_length );
    }

    $bucketizer->optimize(maxrounds => 100_000);

    my @buckets = $bucketizer->buckets();
    $self->count(scalar(@buckets));

    for my $bucket (@buckets) {
        $self->_write_bucket($bucket);
    }

    return 1;
}

sub _write_bucket {
    my $self = shift;
    my $bucket = shift;

    my $staging_dir = $self->temp_staging_directory;
    my $bucket_file = File::Spec->join($staging_dir, $bucket->serial);

    Genome::Sys->write_file($bucket_file, map "$_\n", $bucket->items);
}

1;
