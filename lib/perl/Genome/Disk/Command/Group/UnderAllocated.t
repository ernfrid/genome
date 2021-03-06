#!/usr/bin/env genome-perl

BEGIN {     
    $ENV{UR_DBI_NO_COMMIT} = 1;
}

use strict;
use warnings;

use above "Genome";
use Test::MockObject;
use Test::More;

use_ok('Genome::Disk::Command::Group::UnderAllocated') or die;

subtest 'setup' => sub {
    plan tests => 1;

    my $volume = Test::MockObject->new();
    ok($volume, 'create mock volume');
    $volume->set_always('mount_path', 'obi-wan');
    $volume->set_always('total_kb', 100);
    $volume->mock('percent_used', sub{ 
            my $v = shift;
            return sprintf("%.2f", ( $v->used_kb / $v->total_kb ) * 100);
        });
    $volume->mock('percent_allocated', sub{ 
            my $v = shift;
            return sprintf("%.2f", ( $v->allocated_kb / $v->total_kb ) * 100);
        });

    Sub::Install::reinstall_sub({
            into => 'Genome::Disk::Volume',
            as => 'get',
            code => sub{ return $volume; },
        });

};

subtest 'volume not under allocated' => sub{
    plan tests => 1;

    my $volume = Genome::Disk::Volume->get;
    $volume->set_always('allocated_kb', 50);
    $volume->set_always('used_kb', 5);

    ok( 
        Genome::Disk::Command::Group::UnderAllocated->execute(
            disk_group_names => [ 'jedi' ],
        ),
        'Successfully executed under allocated command',
    );

};

subtest 'volume under allocated' => sub{
    plan tests => 3;

    my $volume = Genome::Disk::Volume->get;
    $volume->set_always('allocated_kb', 50);
    $volume->set_always('used_kb', 55);

    my $cmd = Genome::Disk::Command::Group::UnderAllocated->create(
        disk_group_names => [ 'jedi' ],
    );
    ok($cmd, 'Successfully create under allocated command');
    ok(!$cmd->execute, 'Execute returned undef when disk is under allocated');
    my @msgs = $cmd->status_message;
    is($msgs[0], "Group jedi\n\tVolume obi-wan using 55 kB (55.00 %) but only 50 kB (50.00 %) allocated\n", 'correct status message');

};

done_testing();
