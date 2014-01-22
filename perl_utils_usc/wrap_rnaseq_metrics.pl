#!/usr/bin/perl
use File::Basename;

my $input = $ARGV[0] || die "need input bam file";
my $ref = $ARGV[1] || die "need reference genome file";
my $output = $ARGV[1] || die "need reference genome file";

my $numcores = `cat /proc/cpuinfo | grep processor -c`;
chomp $numcores;
my $ram = 1.3 * $numcores;
my $ram = int($ram);

my $PICARD = "/home/uec-00/shared/production/software/picard/default/CollectRnaSeqMetrics.jar";
my $JAVA = "/home/uec-00/shared/production/software/java/default_java7/bin/java -Xmx$ram" . "G";
my $refDir = "/home/uec-00/shared/production/genomic-data-misc/refSeq_2014_01_21";

## define required file by provided reference genome
$refFlat = "$refDir/refFlathg18.txt" if $ref =~/hg18/ ;
$refFlat = "$refDir/refFlathg19.txt" if $ref =~/hg19/; 
$refFlat = "$refDir/refFlatmm9.txt" if $ref =~/mm9/ ;
$refFlat = "$refDir/refFlatmm10.txt" if $ref =~/mm10/; 

die unless $refFlat;


#my $cmd .= "$JAVA -jar $PICARD REF_FLAT=$refFlat INPUT=$input OUTPUT=$output REFERENCE_SEQUENCE=$ref STRAND_SPECIFICITY=NONE METRIC_ACCUMULATION_LEVEL=ALL_READS TMP_DIR=/export/uec-gs1/laird/shared/tmp";
my $cmd .= "$JAVA -jar $PICARD REF_FLAT=$refFlat INPUT=$input OUTPUT=$output REFERENCE_SEQUENCE=$ref STRAND_SPECIFICITY=NONE TMP_DIR=/export/uec-gs1/laird/shared/tmp";
runcmd($cmd);

sub runcmd{
	my $cmd=shift @_;
	my $caller=(caller(1))[3];
	print STDERR "$caller\t$cmd\n";
	system($cmd);
}
