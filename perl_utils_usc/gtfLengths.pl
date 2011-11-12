#!/usr/bin/env perl

use strict;

foreach my $fn (@ARGV)
{
    my $resa = `wc -l $fn`;
    chomp $resa;
    $resa--;
    my $res = `java GtfFeatLengths $fn 2> /dev/null | tail -1`;
    chomp $res;
    $res =~ s/Total:\s+//g;
    print sprintf("%s:\t%10d\t%6d\n", $fn, $res, $resa);
}
