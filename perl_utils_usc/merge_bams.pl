#!/usr/bin/perl
use strict;
use File::Basename;
use lib dirname (__FILE__);
use EpigenomeUtils;

my $output = shift @ARGV;
my $date = `date`; chomp $date;
my $cmd = "VALIDATION_STRINGENCY=SILENT ASSUME_SORTED=true MERGE_SEQUENCE_DICTIONARIES=true CREATE_INDEX=true USE_THREADING=true MAX_RECORDS_IN_RAM=2000000 OUTPUT='$output' ";

die "merge_bams.pl OUTPUT_FILENAME.bam INPUT1.BAM INPUT2.BAM ...\n" if !@ARGV;

#handle the 1 file case, ie dont merge
if (scalar @ARGV == 1)
{
	die "$ARGV[0] does not exist\n" if !-e $ARGV[0];
	runcmd("cp $ARGV[0] $output");
	runcmd("samtools index $output");
}
else
{
	#check to see if files have readgroups
	my $fileReadGroupCount = 0;
	for my $file (@ARGV)
	{
		$fileReadGroupCount++ if &hasReadGroup($file);
	}

	for my $file (@ARGV)
	{
		die "$file does not exist\n" if !-e $file;
		if($fileReadGroupCount > 0 && $fileReadGroupCount < scalar(@ARGV) && !&hasReadGroup($file))
		{
			#we are catching the case where a non-rg bam is mixed in with a bunch of bams that do have rg. adding the rg to the non-rg bam before merge.
			print STDERR "some files have readgroups and some dont, so adding readgroup to $file since it didnt\n";
			my $RGfile = basename($file) . ".fixed.bam";
			&addReadGroup($file, $RGfile);
			$cmd .= "INPUT='$RGfile' ";
		}
		else
		{
			$cmd .= "INPUT='$file' ";
		}
	}

	runcmd("$JAVA -Xmx4g -jar $PICARD/MergeSamFiles.jar $cmd");
	my $bai = $output;
	$bai =~ s/bam$/bai/;
	runcmd("mv $bai $output.bai");
}

#ResultCount_C0F4UACXX_2_IYE1361A11

&addReadGroup($output, $output) if(! &hasReadGroup($output));

runcmd("$SAMTOOLS flagstat $output > $output\.flagstat\.metric\.txt");

my $outputdups = $output;
$outputdups =~ s/bam$/mdups\.bam/;
runcmd("$JAVA -Xms7g -Xmx7g -jar $PICARD/MarkDuplicates.jar CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT METRICS_FILE=dupmets.txt READ_NAME_REGEX=null INPUT=$output OUTPUT=$outputdups");
my $dupbai = $outputdups;
$dupbai =~ s/bam$/bai/;
runcmd("mv $dupbai $outputdups\.bai");
runcmd("$SAMTOOLS flagstat $outputdups > $outputdups\.flagstat\.metric\.txt");

sub addReadGroup
{
	my $bamIn = shift @_;
	my $bamOut = shift @_;
	my $flowcell;
	my $lane;
	my $lib;
	my $bamInBase = basename($bamIn);;

	if($bamInBase =~ /^(.+?)_(.+?)_(\d+)_(.+?)\./)
	{
		$flowcell = $2;
		$lane = $3;
		$lib = $4;
	}
	elsif($output =~ /^(.+?)_(.+?)_(\d+)_(.+?)\./)
	{
		$flowcell = $2;
		$lane = $3;
		$lib = $4;
	}
	else
	{
		$flowcell = "ANALYSIS";
		$lane = "1";
		$lib = "UNKNOWN_LIB";
	}

	runcmd("$JAVA -Xmx4g -jar $PICARD/AddOrReplaceReadGroups.jar CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT SORT_ORDER=coordinate MAX_RECORDS_IN_RAM=1000000 INPUT='$bamIn' OUTPUT='with_rg_$bamOut' RGID='$flowcell\.$lane' RGLB='$lib' RGPL='illumina Hiseq' RGPU='$flowcell\.$lane' RGSM='$lib' RGCN='USC EPIGENOME CENTER' RGDS='from file $bamIn on $date'");
	#overwrite non-readgroups bams
	runcmd("mv with_rg_$bamOut $bamOut");
	my $bai = "with_rg_$bamOut";
	$bai =~ s/bam$/bai/;
	runcmd("mv $bai $bamOut.bai");
}

sub hasReadGroup
{
	print STDERR "checking for readgroup: ";
	my $bam = shift @_;
	my $header = `$SAMTOOLS view -H $bam`; 
	print "YES, $bam has RG\n" if($header =~ /^\@RG/m);
	print "NO, $bam no RG found\n" if($header !~ /^\@RG/m);
	return 1 if($header =~ /^\@RG/m);
	return 0;
}
