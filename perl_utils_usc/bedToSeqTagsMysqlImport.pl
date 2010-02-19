#!/usr/bin/env perl

use strict;



my $seen = {};
LINE: while (my $line = <STDIN>)
{
    chomp $line;
    my @f = split(/\t/, $line);
    
    # Check valid line while we're getting the strand.
    my $readStart;
    my $valid = 0;
    my $strand = $f[5];
    if ($strand eq "+")
    {
		$valid = 1;
		$readStart = $f[1];
    }
    elsif ($strand eq "-")
    {
		$valid = 1;
		$readStart = $f[2];
    }
    elsif ($strand eq ".")
    {
		$valid = 1;
		$readStart = $f[1];
    }

	my $key = join("__",$readStart,$strand);
	if ($seen->{$key})
	{
		#print STDERR "Already saw $key, skipping ..\n";
		$valid = 0;
	}
	$seen->{$key}++;


    next LINE unless ($valid);

#    print join("\t", $f[0], $readStart, $f[5])."\n";
    print join("\t", $readStart, $strand, 1.0, 1.0, 0, 0, 0, 0, 0,0,0)."\n"; # To fit into methyl schema
}
