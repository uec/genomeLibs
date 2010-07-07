#!/usr/bin/env perl

use strict;
use File::Basename qw/basename/;

my $USAGE = "readDepthFileExtractCounts.pl file1.csv file2.csv ...";


my @fns = @ARGV;

foreach my $fn (@fns)
{
    my $base = basename($fn, qw/.csv/);
    my $outfn = $base."-countcol.csv";

    print STDERR "Working on $fn -> $outfn\n";

    die "Can't read $fn\n" unless (open(IN,$fn));
    die "Can't write to $outfn\n" unless (open(OUT,">$outfn"));

    while (my $line = <IN>)
    {
	my @f = split(/,/,$line);
	my $nC = scalar(@f);
	print OUT $f[$nC-1];
    }

    close(OUT);
    close(IN);
    
}
