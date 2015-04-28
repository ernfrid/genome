#!/usr/bin/env genome-perl

BEGIN { 
    $ENV{UR_DBI_NO_COMMIT} = 1;
    $ENV{UR_USE_DUMMY_AUTOGENERATED_IDS} = 1;
}

use strict;
use warnings;

use above "Genome";
use Test::More;

use_ok('Genome::Disk::Command::Group::AvailableSpace') or die;

my $cmd = Genome::Disk::Command::Group::AvailableSpace->create(
    disk_group_names => Genome::Config::get('disk_group_dev'),
    send_alert => 0,
);
ok($cmd, 'Successfully created avaiable space command object') or die;

my $rv = $cmd->execute;
ok($rv, 'Successfully executed command');

done_testing();
