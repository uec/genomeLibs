#!/usr/bin/env perl

use strict;

while (my $line=<>)
{
    chomp $line;
    next if ($line =~ /^id/i); # header
    my @f = split(/\t/,$line);

    if (scalar(@f)!=6)
    {
	print STDERR "Following line doesn't have 6 cols: ${line}\n";
    }
    else
    {
	my ($id, $chr, $tss, $strand, $s, $e) = @f;

	my $extra = "";
	$extra .= "gene_id \"${id}\" ;transcript_id \"${id}\"; " if ($id);
	$extra .= "mid \"$tss\";";

	print join("\t",
		   $chr,
		   "island",
		   "exon",
		   $s,
		   $e,
		   ".",
		   $strand,
		   ".",
		   $extra) . "\n";
		   
    }
}
