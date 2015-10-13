#!/usr/bin/env perl

my ($chr, $STEP) = @ARGV;

use strict;

my @f = (
    $chr,
    "step$STEP",
    "exon",
    0,
    0,
    ".",
    ".",
    "."
    );

for (my $i = 0; $i <= 250000000; $i+=$STEP)
{
    $f[3] = $i+1;
    $f[4] = $i+$STEP;
    print join("\t",@f)."\n";
}
