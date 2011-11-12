#!/usr/bin/env perl

use strict;
my $FROM = "hg18";
my $TO = "Hg19";

my $libfn = sprintf("/home/rcf-40/bberman/lib/%sTo%s.over.chain",$FROM,$TO);
die "Can't find lib file $libfn\n" unless (-f $libfn);

foreach my $fn (@ARGV)
{

    my $bwc = `wc -l $fn`;
    chomp $bwc;
    $bwc--;


    my $tofn = $fn;
    $tofn =~ s/\.[^\.]*$//g;
    unless ($tofn =~ s/$FROM/$TO/gi)
    {
	$tofn .= ".${TO}";
    }
    $tofn .= ".gtf";
    my $cmd = "liftOver -gff $fn $libfn $tofn unmapped 2> /dev/null ; rm -f unmapped";
#    print STDERR "$cmd\n";
    my $out = `$cmd`;

    my $wc = `wc -l $tofn`;
    chomp $wc;
    $wc--;
    print STDERR sprintf("%s:%7d lines,\t%s:%7d lines\n",$fn, $bwc, $tofn,$wc);
}
