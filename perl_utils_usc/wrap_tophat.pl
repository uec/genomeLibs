#!/usr/bin/perl
$execmd = join(" ", @ARGV);
system($execmd);
$file = scalar(@ARGV) > 3 ? $ARGV[$#ARGV - 1] : $ARGV[$#ARGV];
system("mv coverage.wig " . $file. ".tophat_coverage.wig");
system("mv accepted_hits.sam " . $file. ".tophat_hits.sam");
system("mv junctions.bed " . $file. ".tophat_junctions.bed");


