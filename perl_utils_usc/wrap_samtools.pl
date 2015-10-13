#!/usr/bin/perl
use strict;
use File::Basename;
use lib dirname (__FILE__);
use EpigenomeUtils;

my $samtoolArgs = join ' ', @ARGV;

my $cmd = "$SAMTOOLS $samtoolArgs";
runcmd("$cmd");



sub runcmd
{
        my $cmd = shift @_;
        print STDERR "$cmd\n";
        system($cmd);
}
