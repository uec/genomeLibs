#!/usr/bin/perl
use File::Basename;
use lib dirname (__FILE__);
use EpigenomeUtils;

#CONSTANTS
$igvtools = "$SOFTWAREROOT/igvtools/default/igvtools";


$input = $ARGV[0] || die "specify input file";
$output = $ARGV[1] || die "specify input file";
$ref = $ARGV[2] || die "specify input file";


$genome = "hg19";
$genome = "mm9" if $ref =~ /mm9/;
$genome = "hg18" if $ref =~ /hg18/;
$genome = "rn4" if $ref =~ /rn4/;

$cmd = "$igvtools toTDF  $input $output $genome";
print "$cmd\n";
system($cmd);
