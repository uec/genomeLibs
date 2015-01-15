#!/usr/bin/perl -w
# Mon Jul 12 14:03:06 PDT 2010
use File::Basename;
use lib dirname (__FILE__);
use EpigenomeUtils;

my $usage = "wrap_bwa.pl refFa.fa read1.fq [read2.fq] [outfile name]";
die "$usage\n" unless (@ARGV >= 3);

my $bwa = "$SOFTWAREROOT/bwa/default/bwa";


$refFa = $ARGV[0];
$read1 = $ARGV[1];
my $read1SA = $read1;
$read1SA =~ s/\.fq$//i;
$read1SA .= ".sai";
$outfile = $ARGV[$#ARGV];

die "$outfile already exists" if (-e $outfile);

$outfileSAM = $outfile . ".sam";

# check existence of ref genome
die "reference does not exist.\n" if (! -e $refFa);

# check existence of read sequence files
die "need read sequence files\n" unless ( -e $read1 );

#check phred of read
my $phred = `$SOFTWAREROOT/perl_utils_usc/testFastqQualityScale.pl $read1`;
$phred = $phred =~ /64/ ? "-I" : "";

#aln end 1
my $cmd = join(" ", $bwa, "aln", $phred, "-t 15", $refFa, $read1, "> $read1SA");
runcmd($cmd);

die "need read.sai files\n" unless ( -e $read1SA );


if (@ARGV == 4)
{
	$read2 = $ARGV[2];

	# check existence of read sequence files
	die "need read sequence files\n" unless ( -e $read2 );

	$read2SA = $read2;
	$read2SA =~ s/\.fq$//i;
	$read2SA .= ".sai";

	$cmd = join(" ", $bwa, "aln", $phred, "-t 15", $refFa, $read2, "> $read2SA");
	runcmd($cmd);

	die "need read.sai files\n" unless ( -e $read1SA && -e $read2SA );
	#my $alignCMD = join(" ", $bwa, "sampe -A", $refFa, $read1SA, $read2SA, $read1, $read2, "> $outfileSAM");
	my $alignCMD = join(" ", $bwa, "sampe", $refFa, $read1SA, $read2SA, $read1, $read2, "> $outfileSAM");
	runcmd($alignCMD);
}
else
{
	my $alignCMD = join(" ", $bwa, "samse", $refFa, $read1SA, $read1, "> $outfileSAM");
	runcmd($alignCMD);

}



runcmd("$JAVA -Xmx4g -jar $PICARD/SortSam.jar VALIDATION_STRINGENCY=SILENT SORT_ORDER=coordinate INPUT=$outfileSAM OUTPUT=$outfile\.sorted.bam");
runcmd("$JAVA -Xmx4g -jar $PICARD/ReorderSam.jar VALIDATION_STRINGENCY=SILENT REFERENCE=$refFa INPUT=$outfile\.sorted.bam OUTPUT=$outfile");
