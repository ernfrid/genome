package Genome::InstrumentData::Command::Import::WorkFlow::RetrieveSourcePathFromCgHub;

use strict;
use warnings;

use Genome;

require File::Basename;
require File::Spec;
require Genome::Model::Tools::CgHub::GeneTorrent;
require Genome::Model::Tools::CgHub::Metadata;
require Genome::Model::Tools::CgHub::Query;

class Genome::InstrumentData::Command::Import::WorkFlow::RetrieveSourcePathFromCgHub { 
    is => 'Genome::InstrumentData::Command::Import::WorkFlow::RetrieveSourcePath',
    has => {
        lsf_resource => {
            default_value => Genome::Model::Tools::CgHub::GeneTorrent->__meta__->property_meta_for_name('lsf_resource')->default_value,
        },
    },
    has_output => {
        destination_path => {
            calculate_from => [qw/ working_directory source_path /],
            calculate => q| return File::Spec->join($self->working_directory, File::Basename::basename($source_path).'.bam'); |,
            doc => 'Final destination path.',
        }, 
    },
    has_calculated => {
        metadata_file => {
            calculate_from => [qw/ working_directory /],
            calculate => q| return File::Spec->join($working_directory, 'metadata.xml'); |,
        },
        uuid => {
            calculate_from => [qw/ source_path /],
            calculate => q| return File::Basename::basename($source_path); |,
        },
    },
};

sub _retrieve_path {
    my $self = shift;

    my $retrieve_metadata_path_ok = $self->_retrieve_metadata_path;
    return if not $retrieve_metadata_path_ok;

    my $retrieve_bam_path_ok = $self->_retrieve_bam_path;
    return if not $retrieve_bam_path_ok;

    return 1;
}

sub _retrieve_metadata_path {
    my $self = shift;
    $self->debug_message('Retrieve metadata path from CG Hub...');

    my $uuid = $self->uuid;
    $self->debug_message("UUID: $uuid");
    my $metadata_file = $self->metadata_file;
    $self->debug_message("Metadata file: $metadata_file");

    my $query = Genome::Model::Tools::CgHub::Query->execute(
        uuid => $uuid,
        xml_file => $metadata_file,
    );
    if ( not $query->result ) {
        $self->error_message('Failed to execute cg hub query!');
        return;
    }

    if ( not -s $metadata_file ) {
        $self->error_message("Successfully executed query, but the metadata path does not exist!");
        return;
    }

    $self->debug_message('Retrieve metadata path from CG Hub...done');
    return 1;
}

sub _retrieve_bam_path {
    my $self = shift;
    $self->debug_message('Retrieve source path from CG Hub...');

    my $uuid = $self->uuid;
    $self->debug_message("UUID: $uuid");
    my $destination_path = $self->destination_path;
    $self->debug_message("To: $destination_path");

    my $gene_torrent = Genome::Model::Tools::CgHub::GeneTorrent->execute(
        uuid => $uuid,
        target_path => $destination_path,
    );
    if ( not $gene_torrent->result ) {
        $self->error_message('Failed to execute cg hub gene torrent!');
        return;
    }

    if ( not -s $destination_path ) {
        $self->error_message("Successfully executed gene torrent, but destination path does not exist!");
        return;
    }

    $self->debug_message('Retrieve source path from CG Hub...done');
    return 1;
}

sub _retrieve_source_md5_path {
    my $self = shift;
    $self->debug_message('Retrieve source md5 path from CG Hub...');

    my $metadata = Genome::Model::Tools::CgHub::Metadata->create(
        metadata_file => $self->metadata_file,
    );
    my $bam_file_name = $metadata->bam_file_names;
    my $checksum_type = $metadata->checksum_type_for_file_name($bam_file_name);
    $self->debug_message("Checksum type: $checksum_type");
    if ( not $checksum_type or lc($checksum_type) ne 'md5' ) {
        $self->debug_message('Check sum from metadata is not MD5, skipping storing it.');
        return;
    }

    my $checksum_content = $metadata->checksum_content_for_file_name($bam_file_name);
    $self->debug_message("MD5 from metadata: $checksum_content");
    my $fh = Genome::Sys->open_file_for_writing($self->destination_md5_path);
    $fh->print( join("  ", $checksum_content, $bam_file_name)."\n" );
    $fh->close;

    $self->debug_message('Retrieve source md5 pathfrom CG Hub...done');
    return 1;
}

1;

