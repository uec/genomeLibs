#!/usr/bin/perl
#CONSTANTS
$igvtools = "/home/uec-00/shared/production/software/igvtools/igvtools";


$input = $ARGV[0] || die "specify input file";
$output = $ARGV[1] || die "specify input file";
$genome = $ARGV[2] || die "specify input file";


$cmd = "$igvtools tile -c $input $output $genome";
print "$cmd\n";
system($cmd);
