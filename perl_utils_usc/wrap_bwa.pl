#!/usr/bin/perl -w

use File::Basename;
use lib dirname (__FILE__);
use EpigenomeUtils;

my $usage = "wrap_bwa.pl refFa.fa read1.fq [read2.fq] [outfile name]";
die "$usage\n" unless (@ARGV >= 3);

my $bwa = "$SOFTWAREROOT/bwa/latest/bwa";

$refFa = $ARGV[0];
$read1 = $ARGV[1];
$read2 = $ARGV[2] if @ARGV == 4;
my $read1SA = $read1;
$read1SA =~ s/\.fq$//i;
$read1SA .= ".sai";
$outfile = $ARGV[$#ARGV];
die "$outfile already exists" if (-e $outfile);
$outfileSAM = $outfile . ".sam";

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


if($readLength > 80 && $phred !~ /I/)

{
	print STDERR "doing BWA MEM since readlength is $readLength\n";
	#bwa MEM
	my $cmd = join(" ", $bwa, "mem", "-M -t 27", $refFa, $read1, $read2, "> $outfileSAM");
	runcmd($cmd);
}
else
{
	print STDERR "doing BWA align since readlength is $readLength\n";
	#aln end 1
	my $cmd = join(" ", $bwa, "aln", $phred, "-t 27", $refFa, $read1, "> $read1SA");
	runcmd($cmd);

	die "need read.sai files\n" unless ( -e $read1SA );

	if (@ARGV == 4)
	{
		# check existence of read sequence files
		die "need read sequence files\n" unless ( -e $read2 );

		$read2SA = $read2;
		$read2SA =~ s/\.fq$//i;
		$read2SA .= ".sai";

		$cmd = join(" ", $bwa, "aln", $phred, "-t 27", $refFa, $read2, "> $read2SA");
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
}


runcmd("$JAVA -Xmx4g -jar $PICARD/SortSam.jar VALIDATION_STRINGENCY=SILENT SORT_ORDER=coordinate TMP_DIR=$TMP_DIR INPUT=$outfileSAM OUTPUT=$outfile\.sorted.bam");
runcmd("$JAVA -Xmx4g -jar $PICARD/ReorderSam.jar VALIDATION_STRINGENCY=SILENT TMP_DIR=$TMP_DIR REFERENCE=$refFa INPUT=$outfile\.sorted.bam OUTPUT=$outfile");

