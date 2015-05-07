#!/usr/bin/perl -w

use File::Basename;
use lib dirname (__FILE__);
use EpigenomeUtils;

my $usage = "wrap_bwa.pl refFa.fa read1.fq [read2.fq] [outfile name]";
die "$usage\n" unless (@ARGV >= 3);

$ENV{PATH} .= ":" .  dirname($SAMTOOLS);
$ENV{PATH} .= ":" .  "$SOFTWAREROOT/bwa/latest";

my $bwameth = "$SOFTWAREROOT/bwameth/default/bwameth.py";

$refFa = $ARGV[0];
$read1 = $ARGV[1];
$read2 = $ARGV[2] if @ARGV == 4;
$outfile = $ARGV[$#ARGV];
die "$outfile already exists" if (-e $outfile);
$outfileBWA = $outfile . "bwameth";
$outfileBWABam = $outfile . "bwameth.bam";


#get read length
my $readLength = `cat $read1 | head -n 2 | tail -n 1`;
$readLength = length($readLength);

# check existence of ref genome
die "reference does not exist.\n" if (! -e $refFa);

# check existence of read sequence files
die "need read sequence files\n" unless ( -e $read1 );

#check phred of read
my $phred = `$SOFTWAREROOT/perl_utils_usc/testFastqQualityScale.pl $read1`;
$phred = $phred =~ /64/ ? "-I" : "";


#bwa MEM
my $cmd = "$bwameth --threads 12 --prefix $outfileBWA --reference $refFa $read1 $read2";
runcmd($cmd);


runcmd("$JAVA -Xmx4g -jar $PICARD/SortSam.jar VALIDATION_STRINGENCY=SILENT SORT_ORDER=coordinate TMP_DIR=$TMP_DIR INPUT=$outfileBWABam OUTPUT=$outfile\.sorted.bam");
runcmd("$JAVA -Xmx4g -jar $PICARD/ReorderSam.jar VALIDATION_STRINGENCY=SILENT TMP_DIR=$TMP_DIR REFERENCE=$refFa INPUT=$outfile\.sorted.bam OUTPUT=$outfile");

