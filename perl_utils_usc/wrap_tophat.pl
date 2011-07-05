#!/usr/bin/perl
use File::Basename;
my $SAMTOOLS = "/home/uec-00/shared/production/software/samtools/samtools";
$ENV{'BOWTIE_INDEXES'} = "/home/uec-00/shared/production/genomes/bowtie/";
$ENV{'PATH'} .= ":/home/uec-00/shared/production/software/tophat/default:/home/uec-00/shared/production/software/bowtie/default";
$execmd = "tophat " . join(" ", @ARGV);
print "$execmd\n";
system($execmd);
$file = scalar(@ARGV) > 4 ? $ARGV[$#ARGV - 1] : $ARGV[$#ARGV];
$file = basename($file);
system("mv tophat_out/accepted_hits.bam " . $file. ".tophat_hits.bam");
system("mv tophat_out/junctions.bed " . $file. ".tophat_junctions.bed");
system("mv tophat_out/insertions.bed " . $file. ".tophat_insertions.bed");
system("mv tophat_out/deletions.bed " . $file. ".tophat_deletions.bed");
system("$SAMTOOLS index " . $file. ".tophat_hits.bam");