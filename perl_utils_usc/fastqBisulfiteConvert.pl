#!/usr/bin/perl

use strict;
use File::Basename;
use Getopt::Long;

my $USAGE = "fastqBisulfiteConvert.pl --GtoA seqs1.fastq seqs2.fastq ...";
my $GtoA = 0;
GetOptions ('GtoA' => \$GtoA) || die "$USAGE\n";

my (@files) = @ARGV;
die "$USAGE\n" unless (@files);

foreach my $f (@files)
{
    my $file_num = 1;
    my $lineCount = 0;

    my $outc = fileparse($f);
    my $bsSec = ($GtoA) ? "inSilicoBisulfiteGtoA" : "inSilicoBisulfiteCtoT";
    $outc =~ s/\.(\w+)$/.${bsSec}\.$1/g;
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
	    $line = inSilicoConvert($line, $GtoA);
	    $seq_line = $line;
	}

	$seq_so_far .= $line;

	if ($linecount_within == 4)
	{
	    print OUTC $seq_so_far;
	    
	    $seq_so_far = '';
	    $seq_line = '';

	    $seqs_seen++;
	}
	$lineCount++;
    }

    close (OUTC);    
    close (IN);
}

sub inSilicoConvert
{
    my ($line, $GtoA) = @_;

    if ($GtoA)
    {
	$line =~ s/G/A/g;
	$line =~ s/g/a/g;
    }
    else
    {
	$line =~ s/C/T/g;
	$line =~ s/c/t/g;
    }

    return $line;
}
