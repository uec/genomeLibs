#!/usr/bin/perl

#GLOBALS
$findpeaks = "/home/uec-00/shared/production/software/VancouverShortR/fp4/FindPeaks.jar";

#INPUTS
$inputFile = $ARGV[0] || die "specifiy input BAM file";
$inputFile =~ /^(.+)\.bam$/ || die "specficy BAM file";
$prefix = $1;
$dist = $ARGV[1] || die "specify freg length (ex: 200)";

$cmd = "java -jar $findpeaks -no_warning -name $prefix -input $inputFile -output ./ -aligner sam -dist_type 0 $dist -landerwaterman 0.001 -duplicatefilter  -auto_threshold 0.001";

print "$cmd\n";
