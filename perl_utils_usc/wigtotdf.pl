#!/usr/bin/perl
#CONSTANTS
$igvtools = "/home/uec-00/shared/production/software/igvtools/default/igvtools";


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
