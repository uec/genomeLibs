#!/usr/bin/perl

$input  = $ARGV[0] || die;
$nReads  = $ARGV[1] || die;
$output  = $ARGV[2] || die;

#mult by four since fastq has 4 lines per read.
$nReads = $nReads * 4;

#take for N reads from fastq file
$cmd = "head -n $nReads $input > $output";

print STDERR "$cmd\n";
system $cmd;

