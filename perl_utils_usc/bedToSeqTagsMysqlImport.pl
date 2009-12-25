#!/usr/bin/env perl

use strict;

LINE: while (my $line = <STDIN>)
{
    chomp $line;
    my @f = split(/\t/, $line);
    
    # Check valid line while we're getting the strand.
    my $readStart;
    my $valid = 0;
    if ($f[5] eq "+")
    {
	$valid = 1;
	$readStart = $f[1];
    }
    elsif ($f[5] eq "-")
    {
	$valid = 1;
	$readStart = $f[2];
    }
    elsif ($f[5] eq ".")
    {
	$valid = 1;
	$readStart = $f[1];
    }

    next LINE unless ($valid);

#    print join("\t", $f[0], $readStart, $f[5])."\n";
    print join("\t", $readStart, $f[5], 1.0, 0, 0, 0, 0, 0, 0)."\n"; # To fit into methyl schema
}
