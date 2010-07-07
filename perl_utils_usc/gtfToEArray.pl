#!/usr/bin/env perl

use strict;
use Getopt::Long;

my $USAGE = "gtfToEArray < a.gtf > a.earray.txt";

#GetOptions ('prefix=s' => \$prefix, 'genome=s' => \$genome) || die "$USAGE\n";
#print STDERR "prefix=$prefix\tgenome=$genome\n";

LINE: while (my $line = <>)
{
    next LINE if ($line =~ /^\s*\#/);
    next LINE if ($line =~ /^\s*$/);
    next LINE if ($line =~ /^\s*track/);

    chomp $line;
    my @flds = split(/\t/,$line);
    my ($chr, $src, $feat, $s, $e, $score ,$strand,@rest) = @flds;

    print "${chr}:${s}-${e}\n";
}
