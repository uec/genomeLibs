#!/usr/bin/perl

use strict;
use File::Basename;

my $USAGE = "fastqSplit.pl numberOfFiles seqs1.fastq seqs2.fastq ...";

my ($numberOfFiles, @inputFileList) = @ARGV;

if ( (-f $numberOfFiles) || (@inputFileList == 0) )
{
    die "$USAGE\n";
}

sleep(600);
foreach my $inputFile (@inputFileList)
{
	my $totalFileLength = `wc -l $inputFile`;
	$totalFileLength =~ /^\s*(\d+)\s+.+?$/;
	$totalFileLength = $1;
	chomp $totalFileLength;
	$totalFileLength % 4 == 0 || die "input is not divisible by 4";
	my $totalSeqLength = $totalFileLength / 4;
	my $seqsPerFile = $totalSeqLength / $numberOfFiles;
	open(my $inputHandle, "<$inputFile") || die "can't open input";
	
	for my $i (1..$numberOfFiles)
	{
		my $outputFile = fileparse($inputFile); 
    	$outputFile =~ s/\.(\w+)$/.${i}\.$1/g;
    	open(OUTC, ">$outputFile") || die "Can't write to $outputFile\n";
    	
    	if($i == $numberOfFiles)
    	{
    		while(my $line = <$inputHandle>)
    		{
    			print OUTC $line;
    		}    		
    	}
	    else
	    {
	    	for my $j (1..($seqsPerFile * 4))
	    	{
	    		my $line = <$inputHandle>;
	    		print OUTC $line;
	    	}
    	}
    	close(OUTC);    	    	
	}
	close($inputHandle);	
}
