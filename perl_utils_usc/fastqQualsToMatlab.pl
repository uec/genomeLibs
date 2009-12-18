#!/usr/bin/perl

use strict;
use File::Basename;
use Getopt::Long;

my $USAGE = "fastqQualsToMatlab.pl --illuminaQuals seqs1.fastq seqs2.fastq ...";
my $illuminaQuals = 0;
GetOptions ('illuminaQuals' => \$illuminaQuals) || die "$USAGE\n";

my (@files) = @ARGV;
die "$USAGE\n" unless (@files);

foreach my $f (@files)
{
    my $file_num = 1;
    my $lineCount = 0;

    my $outc = fileparse($f);
    my $matlabSec = "matlab-" . ($illuminaQuals?"illuminaQuals":"sangerQuals");
    $outc =~ s/\.(\w+)$/.${matlabSec}\.$1/g;
    die "Can't write to $outc\n" unless (open(OUTC,">$outc"));

    die "Can't read $f\n" unless (open(IN,$f));

    my $seqs_seen = 0;

    my $linecount_within = 0;
    my $seq_so_far = '';
    my $seq_line = '';
    while (my $line = <IN>)
    {
	if ($line =~ /^\@/ && ($lineCount % 4 == 0))
	{
	    $linecount_within = 1;
	}
	else
	{
	    $linecount_within++;
	}
	
	if ($linecount_within == 2)
	{
		my $qualsLine = getQualsLine($line,$illuminaQuals);
		print $outc $qualsLine . "\n";
	}

	$seq_so_far .= $line;

	if ($linecount_within == 4)
	{
	    
	    $seq_so_far = '';
	    $seq_line = '';

	    $seqs_seen++;
	}
	$lineCount++;
    }

    close (OUTC);    
    close (IN);
}

sub getQualsLine
{
    my ($line, $illuminaQuals) = @_;
	my @flds = ();
 

    return join(",",@flds);
}
