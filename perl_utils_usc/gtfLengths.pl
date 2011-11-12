#!/usr/bin/env perl

use strict;

foreach my $fn (@ARGV)
{
    my $resa = `wc -l $fn`;
    chomp $resa;
    $resa--;

    # Feat lengths
    my $res = `java GtfFeatLengths $fn 2> /dev/null | tail -1`;
    chomp $res;
    $res =~ s/Total:\s+//g;

    # Percentage of human genome
    my $perc = 100* ($res / 2.65E9);


    print sprintf("%s:\t%10d\t%7.4f%%\t%6d\n", $fn, $res, $perc, $resa);
}
