#!/usr/bin/perl
use strict;

my $SAMTOOLS = "/home/uec-00/shared/production/software/samtools/samtools";
my $samtoolArgs = join ' ', @ARGV;

my $cmd = "$SAMTOOLS $samtoolArgs";
runcmd("$cmd");



sub runcmd
{
        my $cmd = shift @_;
        print STDERR "$cmd\n";
        system($cmd);
}
