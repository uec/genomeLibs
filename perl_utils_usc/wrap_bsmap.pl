#!/usr/bin/perl
use strict;
use File::Basename;
# bsmap -a ./s_7_1_sequence.200k.txt -d ~/genomes/hg18_unmasked/hg18_unmasked.plusContam.fa -o s_7_1_sequence.200k.sam -z @ -p 11 -b s_7_2_sequence.200k.txt -s 18 -v 10 -q 2

my $BSMAP = "/home/uec-00/shared/production/software/bsmap/default/bsmap";

#change to AMD version in case SSE4.2 not supported on this machine.
$BSMAP = "/home/uec-00/shared/production/software/bsmap/default-amd/bsmap" if `bash -c \"$BSMAP 2>&1\"` =~ /Fatal/s;
print STDERR "$BSMAP\n";
my $SAMTOOLS = "/home/uec-00/shared/production/software/samtools/samtools";
my $PICARD = "/home/uec-00/shared/production/software/picard/default/";
my $JAVA = "/home/uec-00/shared/production/software/java/default/bin/java";
my $MAQ = "/home/uec-00/shared/production/software/maq-0.7.1/maq";

my $output = $ARGV[0] || die "no output specified";
my $genome = $ARGV[1] || die "no genome specified";
my $read1 = $ARGV[2] || die "no read1 specified";
my $read2 = $ARGV[3]; 
my $samOutput = "$output" . ".sam";

#check phred of read

my $read1sanger = $read1; 
my $read2sanger = $read2 if $read2;

#handle compressed inputs
runcmd("zcat $read1 >" . basename($read1) . ".nogz.fastq") if($read1 =~ /gz$/);
runcmd("zcat $read2 >" . basename($read2) . ".nogz.fastq") if($read2 =~ /gz$/);
$read1 = basename($read1) . ".nogz.fastq" if($read1 =~ /gz$/);
$read2 = basename($read2) . ".nogz.fastq" if($read2 =~ /gz$/);

my $phred = `/home/uec-00/shared/production/software/perl_utils_usc/testFastqQualityScale.pl $read1`;
if($phred =~ /64/)
{
	
	#hack, not producing correct phreds in bam
	$read1sanger = basename($read1) . ".fastq";
	$read2sanger = basename($read2) . ".fastq" if $read2;
	runcmd("$MAQ ill2sanger $read1 $read1sanger");
	runcmd("$MAQ ill2sanger $read2 $read2sanger") if $read2;
}
else
{
	$read1sanger = $read1;
	$read2sanger = $read2;
}

my $cmd = "$BSMAP -a $read1sanger ";
$cmd .= "-b $read2sanger " if $read2;
$cmd .= "-d $genome -o $samOutput ";
$cmd .= "-p 12 -s 16 -v 10 -q 2 -A AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT";


runcmd($cmd);

#samtobamsort
runcmd("$JAVA -Xmx7g -jar $PICARD/SortSam.jar VALIDATION_STRINGENCY=SILENT INPUT=$samOutput OUTPUT=$output SORT_ORDER=coordinate TMP_DIR=/export/uec-gs1/laird/shared/tmp");

#samtools sort
#runcmd("$SAMTOOLS view -b -S -o $output $samOutput");


#sort
#print STDERR "sorting bam using samtools\n";
#system("$SAMDIR/samtools sort $output.presort.bam $output.sorted");
#system("mv $output.sorted.bam $output");

#mark duplicates
#my $alldups = $output . ".dups.bam";
#my $alldupsmetrics = $output . ".dups.bam.metrics";
#runcmd("$JAVA -Xmx14g -jar $PICARD/MarkDuplicates.jar VALIDATION_STRINGENCY=SILENT MAX_RECORDS_IN_RAM=3000000 INPUT='$output' OUTPUT='$alldups' METRICS_FILE='$alldupsmetrics'");
#system("mv $alldups $output");

sub runcmd
{
        my $cmd = shift @_;
        print STDERR "$cmd\n";
        system($cmd);
}
