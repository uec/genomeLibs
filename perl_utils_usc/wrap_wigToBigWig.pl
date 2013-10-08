#!/usr/bin/perl

$ucsc = "/auto/uec-00/shared/production/software/UCSC_Browser_Tools/default";

$input = $ARGV[0] || die "specifiy wig file";
$input =~ /wig/i || die "specifiy wig file";
$output = $ARGV[1] || die "specifiy output file";

$genome = $ARGV[2] || "hg19";
$genome = "mm9" if $input =~ /mm9/;
$genome = "mm10" if $input =~ /mm10/;
$genome = "hg18" if $input =~ /hg18/;
$genome = "rn5" if $input =~ /rn5/;
$genome = "danRer7" if $input =~ /zv9/i;

system("$ucsc/wigToBigWig -clip $input $ucsc/$genome\.chrom.sizes $output");
