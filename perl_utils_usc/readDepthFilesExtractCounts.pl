#!/usr/bin/env perl

use strict;
use File::Basename qw/basename/;

my $USAGE = "readDepthFileExtractCounts.pl file1.csv file2.csv ... > out.csv";


my @fns = @ARGV;

my $filenum = 1;
my $lanes = [];
foreach my $fn (@fns)
{
    my $base = basename($fn, qw/.csv/);

    print STDERR "Working on $fn\n";

    die "Can't read $fn\n" unless (open(IN,$fn));

    my @outlane = ();
    my @outcoords = ();
    while (my $line = <IN>)
    {
	chomp $line;

	my @f = split(/,/,$line);
	my $secondStrand = ($f[5] == 1); # "1" always comes second


	if (($filenum==1) && $secondStrand)
	{
	    my $coord= $f[6];
	    push(@outcoords,$coord);
	}
    
	my $nC = scalar(@f);
	my $count = $f[$nC-1];
	if ($secondStrand)
	{
	    @outlane[scalar(@outlane)-1]+=$count;
	}
	else
	{
	    push(@outlane,$count);
	}
    }

    @{$lanes}[0] = \@outcoords if (@outcoords);
    @{$lanes}[$filenum] = \@outlane;

    close(IN);
    $filenum++;
}

my $nL = scalar(@$lanes);
my $nR = scalar(@{@{$lanes}[1]});
print STDERR "nL=$nL, nR=$nR\n";

for (my $i = 0; $i < $nR; $i++)
{
    for (my $j = 0; $j < $nL; $j++)
    {
	my $v = @{@{$lanes}[$j]}[$i];
	print $v;
	print "," unless ($j==($nL-1));
    }
    print "\n";
}
