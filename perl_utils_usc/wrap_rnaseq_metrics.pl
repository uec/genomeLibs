#!/usr/bin/perl
use strict;
use File::Basename;
use lib dirname (__FILE__);
use EpigenomeUtils;

my $input = $ARGV[0] || die "need input bam file";
my $ref = $ARGV[1] || die "need reference genome file";
my $output = $ARGV[2] || die "need reference genome file";

my $numcores = `cat /proc/cpuinfo | grep processor -c`;
chomp $numcores;
my $ram = 1.3 * $numcores;
my $ram = int($ram);


my $refDir = "$GENOMEROOT/genomic-data-misc/refSeq_2014_01_21";

## define required file by provided reference genome
my $refFlat;
$refFlat = "$refDir/refFlathg18.txt" if $ref =~/hg18/ ;
$refFlat = "$refDir/refFlathg19.txt" if $ref =~/hg19/; 
$refFlat = "$refDir/refFlatmm9.txt" if $ref =~/mm9/ ;
$refFlat = "$refDir/refFlatmm10.txt" if $ref =~/mm10/; 

die unless $refFlat;

my $cmd .= "$JAVA -jar $PICARD REF_FLAT=$refFlat INPUT=$input OUTPUT=$output REFERENCE_SEQUENCE=$ref STRAND_SPECIFICITY=NONE VALIDATION_STRINGENCY=SILENT TMP_DIR=$TMP_DIR";
runcmd($cmd);

