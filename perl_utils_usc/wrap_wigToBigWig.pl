#!/usr/bin/perl
use strict;
use File::Basename;
use lib dirname (__FILE__);
use EpigenomeUtils;

my $ucsc = "$SOFTWAREROOT/UCSC_Browser_Tools/default";

my $input = $ARGV[0] || die "specifiy wig file";
$input =~ /wig/i || die "specifiy wig file";
my $output = $ARGV[1] || die "specifiy output file";

my $genome = $ARGV[2] || "hg19";
$genome = "mm9" if $input =~ /mm9/;
$genome = "mm10" if $input =~ /mm10/;
$genome = "hg18" if $input =~ /hg18/;
$genome = "rn5" if $input =~ /rn5/;
$genome = "danRer7" if $input =~ /zv9/i;


system("$ucsc/wigToBigWig -clip $input $ucsc/$genome\.chrom.sizes $output");
