#!/usr/bin/perl
use File::Basename;
my $SAMTOOLS = "/home/uec-00/shared/production/software/samtools/samtools";
#$ENV{'BOWTIE_INDEXES'} = "/home/uec-00/shared/production/genomes/bowtie/";
$ENV{'PATH'} .= ":/home/uec-00/shared/production/software/tophat2/default:/home/uec-00/shared/production/software/bowtie2/default";

$file = scalar(@ARGV) > 4 ? $ARGV[$#ARGV - 1] : $ARGV[$#ARGV];
$file = basename($file);


my $phred = `/home/uec-00/shared/production/software/perl_utils_usc/testFastqQualityScale.pl $file`;
$phred = $phred =~ /64/ ? "--solexa1.3-quals" : "";


$execmd = "tophat2 $phred --no-coverage-search " . join(" ", @ARGV);
print "$execmd\n";
system($execmd);
system("mv tophat_out/accepted_hits.bam " . $file. ".tophat_hits.bam");
system("mv tophat_out/unmapped.bam " . $file. ".unmapped.bam");
system("mv tophat_out/junctions.bed " . $file. ".tophat_junctions.bed");
system("mv tophat_out/insertions.bed " . $file. ".tophat_insertions.bed");
system("mv tophat_out/deletions.bed " . $file. ".tophat_deletions.bed");
system("$SAMTOOLS merge tmp.bam " . $file. ".tophat_hits.bam " . $file. ".unmapped.bam");
system("mv tmp.bam " . $file. ".tophat_hits.bam");
system("$SAMTOOLS index " . $file. ".tophat_hits.bam");
