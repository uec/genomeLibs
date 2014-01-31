#!/usr/bin/perl

use strict;
use File::Basename;

my $USAGE = "fastqSplit.pl numberOfFiles seqs1.fastq seqs2.fastq ...";
my $seqret = "/home/uec-00/shared/production/software/emboss/default/bin/seqret";

my ($numberOfFiles, @inputFileList) = @ARGV;

open(METRIC, ">inputreads.metric.txt");

if ( (-f $numberOfFiles) || (@inputFileList == 0) )
{
    die "$USAGE\n";
}

if ($numberOfFiles == 1)
{
	foreach my $inputFile (@inputFileList)
	{
		 if($inputFile =~ /bz2$/)
		 {
		 	my $unzippedInputFile = basename($inputFile);
		 	$unzippedInputFile =~ s/\.bz2$//;
			system("bzip2 -d -k -c $inputFile > $unzippedInputFile") ;
			$inputFile = $unzippedInputFile;
		 }
		 if($inputFile =~ /gz$/)
		 {
		 	my $unzippedInputFile = basename($inputFile);
		 	$unzippedInputFile =~ s/\.gz$//;
			my $cmd = "gzip -d -c $inputFile > $unzippedInputFile";
			print STDERR $cmd;
			system($cmd) ;
			#print("$seqret fastq-sanger::$unzippedInputFile fastq-illumina:$unzippedInputFile\.tmp");
			#system("$seqret fastq-sanger::$unzippedInputFile fastq-illumina:$unzippedInputFile\.tmp");
			#system("mv $unzippedInputFile\.tmp $unzippedInputFile");
			$inputFile = $unzippedInputFile;
		 }
		 my $outputFile = basename($inputFile);
		 $outputFile =~ s/\.(\w+)$/\.1\.$1/g;
		 print STDERR "cp $inputFile $outputFile\n";
		 system("cp $inputFile $outputFile");
	}
	exit;
}


foreach my $inputFile (@inputFileList)
{
	if($inputFile =~ /bz2$/)
	 {
		my $unzippedInputFile = basename($inputFile);
		$unzippedInputFile =~ s/\.bz2$//;
		system("bzip2 -d -k -c $inputFile > $unzippedInputFile") ;
		$inputFile = $unzippedInputFile;
	 }
	my $totalFileLength = `wc -l $inputFile`;
	$totalFileLength =~ /^\s*(\d+)\s+.+?$/;
	$totalFileLength = $1;
	chomp $totalFileLength;
	$totalFileLength % 4 == 0 || die "input is not divisible by 4";
	my $totalSeqLength = $totalFileLength / 4;
	my $seqsPerFile = $totalSeqLength / $numberOfFiles;
	$seqsPerFile =~ s/\.\d+$//;
	print "length:  $totalFileLength\n";
	print "spf:  $seqsPerFile\n";
	print METRIC "$inputFile\t$totalSeqLength\n";

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
