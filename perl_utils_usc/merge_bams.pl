#!/usr/bin/perl

my $SAMDIR = "/home/uec-00/shared/production/software/samtools";
my $PICARD = "/home/uec-00/shared/production/software/picard/default/";
my $JAVA = "/home/uec-00/shared/production/software/java/default/bin/java";

my $output = shift @ARGV;

my $cmd = "VALIDATION_STRINGENCY=SILENT MERGE_SEQUENCE_DICTIONARIES=true CREATE_INDEX=true USE_THREADING=true MAX_RECORDS_IN_RAM=3000000 OUTPUT='$output' ";
for my $file (@ARGV)
{
	$cmd .= "INPUT='$file' ";
}

runcmd("$JAVA -Xmx14g -jar $PICARD/MergeSamFiles.jar $cmd");
sub runcmd
{
        my $cmd = shift @_;
        print STDERR "$cmd\n";
        system($cmd);
}
