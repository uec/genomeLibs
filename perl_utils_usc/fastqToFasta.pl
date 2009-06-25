#!/usr/bin/perl

use strict;

my $USAGE = "fastqToFasta.pl seqs1.fastq seqs2.fastq ...\nfastqToFasta.pl < seqs1.fastq > seqs1.fa\n";

my (@files) = @ARGV;

if (@files == 0)
{
    process(\*STDIN, \*STDOUT);
}
else
{
    foreach my $f (@files)
    {
	my $outc = $f; $outc =~ s/\.(\w+)$/.fa/g;
	die "Can't write to $outc\n" unless (open(OUTC,">$outc"));
	die "Can't read $f\n" unless (open(IN,$f));
	process(\*IN, \*OUTC);
	close (OUTC);    
	close (IN);
    }
}


sub process
{
    my ($inFh, $outFh) = @_;

    my $seqs_seen = 0;

    my $linecount_within = 0;
    my $seq_so_far = '';
    my $seq_line = '';
    my $seq_name = '';
    while (my $line = <$inFh>)
    {
	if ($line =~ /^\@(.*)$/)
	{
	    $linecount_within = 1;
	    $seq_name = $1;
	}
	else
	{
	    $linecount_within++;
	}
	
	$seq_so_far .= $line;
	
	if ($linecount_within == 2)
	{
	    $seq_line = $line;
	}
	elsif ($linecount_within == 4)
	{
	    print $outFh ">${seq_name}\n$seq_line";
	    
	    $seq_so_far = '';
	    $seq_line = '';
	    $seq_name = '';

	    $seqs_seen++;
	}
    }
}
