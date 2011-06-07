#!/usr/bin/perl

my $SAMTOOLS = "/home/uec-00/shared/production/software/samtools/samtools";
my $PICARD = "/home/uec-00/shared/production/software/picard/default/";
my $JAVA = "/home/uec-00/shared/production/software/java/default/bin/java";

my $output = shift @ARGV;

my $cmd = "VALIDATION_STRINGENCY=SILENT ASSUME_SORTED=true  MERGE_SEQUENCE_DICTIONARIES=true CREATE_INDEX=true USE_THREADING=true MAX_RECORDS_IN_RAM=3000000 OUTPUT='$output' ";
for my $file (@ARGV)
{
        $cmd .= "INPUT='$file' ";
}
runcmd("$JAVA -Xmx14g -jar $PICARD/MergeSamFiles.jar $cmd");
my $bai = $output;
$bai =~ s/bam$/bai/;
runcmd("mv $bai $output.bai");
runcmd("$SAMTOOLS flagstat $output > $output\.flagstat\.metric\.txt");

$outputdups = $output;
$outputdups =~ s/bam$/mdups\.bam/;
runcmd("$JAVA -Xms14g -Xmx14g -jar $PICARD/MarkDuplicates.jar CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT METRICS_FILE=dupmets.txt READ_NAME_REGEX=null INPUT=$output OUTPUT=$outputdups");
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