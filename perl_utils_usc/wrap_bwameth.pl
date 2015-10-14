#!/usr/bin/perl -w

use File::Basename;
use lib dirname (__FILE__);
use EpigenomeUtils;

my $usage = "wrap_bwa.pl refFa.fa read1.fq [read2.fq] [outfile name]";
die "$usage\n" unless (@ARGV >= 3);

$ENV{PATH} .= ":" .  dirname($SAMTOOLS);
$ENV{PATH} .= ":" .  "$SOFTWAREROOT/bwa/latest";
$ENV{PATH} = "$SOFTWAREROOT/python/default/bin:" . $ENV{PATH};

my $bwameth = "$SOFTWAREROOT/bwameth/default/bwameth.py";
my $date = `date`; chomp $date;

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
my $cmd = "$bwameth --threads 27 --prefix $outfileBWA --reference $refFa $read1 $read2";
runcmd($cmd);


runcmd("$JAVA -Xmx4g -jar $PICARD/SortSam.jar VALIDATION_STRINGENCY=SILENT SORT_ORDER=coordinate TMP_DIR=$TMP_DIR INPUT=$outfileBWABam OUTPUT=$outfile\.sorted.bam");
my $flowcell = "flowcell$date";
my $lane = "1";
my $lib = basename($outfile);
$lib =~ s/\W//g;
runcmd("$JAVA -Xmx4g -jar $PICARD/AddOrReplaceReadGroups.jar CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT SORT_ORDER=coordinate MAX_RECORDS_IN_RAM=1000000 INPUT='$outfile\.sorted.bam' OUTPUT='$outfile\.withrg.bam' RGID='$flowcell\.$lane' RGLB='$lib' RGPL='illumina' RGPU='$flowcell\.$lane' RGSM='$lib' RGCN='VARI' RGDS='from file $outfile on $date' TMP_DIR=$TMP_DIR");
runcmd("$JAVA -Xmx4g -jar $PICARD/ReorderSam.jar VALIDATION_STRINGENCY=SILENT TMP_DIR=$TMP_DIR REFERENCE=$refFa INPUT=$outfile\.withrg.bam OUTPUT=$outfile");

