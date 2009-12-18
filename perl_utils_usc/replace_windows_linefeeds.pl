#!/usr/bin/perl

use strict;

while (my $line = <>)
{
    chomp;
    my @lines_fixed = split(chr(0xd), $line);
    print join("\n",grep {/\S+/} @lines_fixed);
    print "\n";
}
