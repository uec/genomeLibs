#!/usr/bin/perl

use strict;

my $USAGE = "fastqSimpleRepeatAnalysis.pl GAATG GAATGGAATG CATTC CATTCCATTC file1.fastq file2.fastq";

my $startedFiles = 0;
my @patterns = ();
my @files = ();
while (my $arg = shift(@ARGV))
{
    if (-f $arg)
    {
	push(@files, $arg);
	$startedFiles = 1;
    }
    else
    {
	die "$USAGE\n" if ($startedFiles);
	push(@patterns, $arg);
    }
}

die "$USAGE\n" unless (@files && @patterns);


my $nPats = scalar(@patterns);

foreach my $f (@files)
{
    # Open input
    die "Can't read $f\n" unless (open(IN,$f));

    my $linecount_within = 1;
    my $linecount_global = 1;
    my $seqs_seen = 0;
    my $seq_so_far = '';
    my $seq_line = '';
    my $seqLen = 0;
    my $arrays = [];
    while (my $line = <IN>)
    {
	$seq_so_far .= $line;
	
	if ($linecount_within == 1)
	{
	    # Double check the format
	    print STDERR "Incorrect FASTQ file $f\nLine ${linecount_global}: $line\nMod4 lines should start with \@\n"
		unless ($line =~ /^\@/);
	}
	elsif ($linecount_within == 2)
	{
	    $seq_line = $line;
	}
	elsif ($linecount_within == 4)
	{
	    chomp $seq_line;

	    # Initialize arrays if it's the first seq
	    if (!$seqs_seen)
	    {
		$seqLen = length($seq_line);
		for (my $i=0; $i<$nPats; $i++)
		{
		    @{$arrays}[$i] = [map {0} (1..$seqLen)];
		}
	    }
	    
	    # Increment arrays
	    for (my $i=0; $i<$nPats; $i++)
	    {
		
	    }
	    
	    # Increment other
	    $seqs_seen++;

	    # Reset our linecount
	    $seq_so_far = '';
	    $seq_line = '';
	    $linecount_within = 0;
	}

	# Increment
	$linecount_within++;
	$linecount_global++;
    }

    # Output the hashes.
    
    close (IN);
}


sub analyzeSequence
{
    my ($fn, $pats, $seq) = @_;

    

}
