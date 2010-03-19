#!/usr/bin/perl
$ENV{'BOWTIE_INDEXES'} = "/home/uec-00/shared/production/genomes/bowtie/";
$execmd = join(" ", @ARGV);
system($execmd);
$file = scalar(@ARGV) > 3 ? $ARGV[$#ARGV - 1] : $ARGV[$#ARGV];
system("mv tophat_out/coverage.wig " . $file. ".tophat_coverage.wig");
system("mv tophat_out/accepted_hits.sam " . $file. ".tophat_hits.sam");
system("mv tophat_out/junctions.bed " . $file. ".tophat_junctions.bed");


