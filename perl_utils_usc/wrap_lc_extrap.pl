#!/usr/bin/perl
use File::Basename;
use lib dirname (__FILE__);
use EpigenomeUtils;


#GLOBALS
$lc_extrap = "$SOFTWAREROOT/AndrewSmithTools/preseq-1.0.2.Linux_x86_64/preseq lc_extrap";

#INPUTS
$inputFile = $ARGV[0] || die "specifiy input BAM file";
$outputFile = $ARGV[1] || die "specifiy output TXT file";
$inputFile =~ /^(.+)\.bam$/ || die "specficy outout";

$cmd = "$lc_extrap -v -s 10000000 -e 3000000000 -bam $inputFile > $outputFile ";
print "$cmd\n";
system("$cmd");

