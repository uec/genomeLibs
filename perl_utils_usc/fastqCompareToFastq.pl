#!/usr/bin/perl

use strict;
use File::Basename;
use Getopt::Long;

my $USAGE = "fastqCompareToFastq.pl --bisulfiteConvert --cycleStart 0 --cycleEnd 10  seqs1.fastq seqs2.fastq";

my $bs = 0;
my $cycStart = 0;
my $cycEnd = 0;
GetOptions ('bisulfiteConvert' => \$bs, 'cycleStart=i' => \$cycStart, 'cycleEnd=i' => \$cycEnd) || die "$USAGE\n";

die "$USAGE\n" unless (@ARGV==2);
my ($filea, $fileb) = @ARGV;


    my $file_num = 1;
    my $lineCount = 0;

    die "Can't read $filea\n" unless (open(FA,$filea));
    die "Can't read $fileb\n" unless (open(FB,$fileb));

    my $seqs_seen = 0;

    my $linecount_within = 0;
    while ( ($seqs_seen < 1E6) && (my $line = <FA>) )
    {
	my $lineb = <FB>;
	chomp $line;
	chomp $lineb;

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
	    #$line = inSilicoConvert($line, 0);
	    #$lineb = inSilicoConvert($lineb, 1);

	    my $matches = countMatches($line,$lineb, $bs, $cycStart, $cycEnd);
	    
#	    print "$matches\t$line\t$lineb\n";
	    print "$matches\n";

	    print STDERR "On seq $seqs_seen...\n" if (!($seqs_seen % 1E5));
	}

	if ($linecount_within == 4)
	{
	    $seqs_seen++;
	}
	$lineCount++;
    }

    close (FA);
    close (FB);


sub countMatches
{
    my ($seqa, $seqb, $bs, $cycStart, $cycEnd) = @_;

    my $match = 0;
    my $total = 0;
    $cycStart = 0 unless ($cycStart);
    $cycEnd = length($seqa) unless ($cycEnd);
    for (my $i = $cycStart; $i < $cycEnd; $i++)
    {
	my $chara = uc(substr($seqa, $i,1));
	my $charb = uc(substr($seqb, $i,1));
	
#	print "\t$chara\t$charb\n";
	my $isN = ($chara eq 'N') || ($charb eq 'N');
	if (!$isN)
	{	    
	    $match++ if (($chara eq $charb) || ($bs && ($chara eq 'T') && ($charb eq 'C')) ||
			 ($bs && ($chara eq 'G') && ($charb eq 'A')) );
	    $total++;
	}
    }

    my $frac = ($total) ? ($match/$total) : 0;
    return sprintf("%0.2f",$frac);
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
