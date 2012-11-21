#!/usr/bin/perl

my $SAMTOOLS = "/home/uec-00/shared/production/software/samtools/samtools";
my $PICARD = "/home/uec-00/shared/production/software/picard/default/";
my $JAVA = "/home/uec-00/shared/production/software/java/default/bin/java";

my $output = shift @ARGV;
my $date = `date`; chomp $date;

my $cmd = "VALIDATION_STRINGENCY=SILENT ASSUME_SORTED=true MERGE_SEQUENCE_DICTIONARIES=true CREATE_INDEX=true USE_THREADING=true MAX_RECORDS_IN_RAM=2000000 OUTPUT='$output' ";

if (scalar @ARGV == 1)
{
	runcmd("cp $ARGV[0] $output");
	runcmd("samtools index $output");
}
else
{
	for my $file (@ARGV)
	{
		$cmd .= "INPUT='$file' ";
	}
	runcmd("$JAVA -Xmx12g -jar $PICARD/MergeSamFiles.jar $cmd");
	my $bai = $output;
	$bai =~ s/bam$/bai/;
	runcmd("mv $bai $output.bai");
}

#ResultCount_C0F4UACXX_2_IYE1361A11

my $header = `$SAMTOOLS view -H $output`; 
if($output =~ /^(.+?)_(.+?)_(\d+)_(.+?)\./ && $header !~ /^@RG/)
{
	my $flowcell = $2;
	my $lane = $3;
	my $lib = $4;
	runcmd("$JAVA -Xmx4g -jar $PICARD/AddOrReplaceReadGroups.jar CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT SORT_ORDER=coordinate MAX_RECORDS_IN_RAM=1000000 INPUT='$output' OUTPUT='with_rg_$output' RGID='$flowcell\.$lane' RGLB='$lib' RGPL='illumina Hiseq' RGPU='$flowcell\.$lane' RGSM='$lib' RGCN='USC EPIGENOME CENTER' RGDS='from file $output on $date'");
	#overwrite non-readgroups bams
	runcmd("mv with_rg_$output $output");
	my $bai = "with_rg_$output";
	$bai =~ s/bam$/bai/;
	runcmd("mv $bai $output.bai");

}



runcmd("$SAMTOOLS flagstat $output > $output\.flagstat\.metric\.txt");

$outputdups = $output;
$outputdups =~ s/bam$/mdups\.bam/;
runcmd("$JAVA -Xms12g -Xmx12g -jar $PICARD/MarkDuplicates.jar CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT METRICS_FILE=dupmets.txt READ_NAME_REGEX=null INPUT=$output OUTPUT=$outputdups");
my $dupbai = $outputdups;
$dupbai =~ s/bam$/bai/;
runcmd("mv $dupbai $outputdups\.bai");
runcmd("$SAMTOOLS flagstat $outputdups > $outputdups\.flagstat\.metric\.txt");


sub runcmd
{
        my $cmd = shift @_;
        print STDERR "$cmd\n";
        system($cmd);
}
