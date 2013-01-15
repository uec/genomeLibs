#!/usr/bin/perl

#GLOBALS
$lc_extrap = "/home/uec-00/shared/production/software/AndrewSmithTools/RationalFunctionComplexity/lc_extrap";

#INPUTS
$inputFile = $ARGV[0] || die "specifiy input BAM file";
$outputFile = $ARGV[1] || die "specifiy output TXT file";
$inputFile =~ /^(.+)\.bam$/ || die "specficy outout";

$cmd = "$lc_extrap -v -s 10000000 -bam $inputFile > $outputFile ";
print "$cmd\n";
system("$cmd");

