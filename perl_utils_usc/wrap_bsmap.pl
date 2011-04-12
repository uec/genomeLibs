#!/usr/bin/perl
use strict;
# bsmap -a ./s_7_1_sequence.200k.txt -d ~/genomes/hg18_unmasked/hg18_unmasked.plusContam.fa -o s_7_1_sequence.200k.sam -z @ -p 11 -b s_7_2_sequence.200k.txt -s 18 -v 10 -q 2

my $BSMAP = "/home/uec-00/shared/production/software/bsmap/default/bsmap";
my $SAMTOOLS = "/home/uec-00/shared/production/software/samtools/samtools";
my $PICARD = "/home/uec-00/shared/production/software/picard/default/";
my $JAVA = "/home/uec-00/shared/production/software/java/default/bin/java";

my $output = $ARGV[0] || die "no output specified";
my $genome = $ARGV[1] || die "no genome specified";
my $read1 = $ARGV[2] || die "no read1 specified";
my $read2 = $ARGV[3]; 

my $samOutput = "$output" . ".sam";


my $cmd = "$BSMAP -a $read1 ";
$cmd .= "-b $read2 " if $read2;
$cmd .= "-d $genome -o $samOutput ";
$cmd .= "-z \@ -p 8 -s 16 -v 10 -q 2";

#runcmd($cmd);

#samtobamsort
runcmd("$JAVA -Xmx14g -jar $PICARD/SortSam.jar VALIDATION_STRINGENCY=SILENT INPUT=$samOutput OUTPUT=$output SORT_ORDER=coordinate");

#samtools sort
#runcmd("$SAMTOOLS view -b -S -o $output $samOutput");


#sort
#print STDERR "sorting bam using samtools\n";
#system("$SAMDIR/samtools sort $output.presort.bam $output.sorted");
#system("mv $output.sorted.bam $output");

#mark duplicates
my $alldups = $output . ".dups.bam";
my $alldupsmetrics = $output . ".dups.bam.metrics";
runcmd("$JAVA -Xmx14g -jar $PICARD/MarkDuplicates.jar VALIDATION_STRINGENCY=SILENT MAX_RECORDS_IN_RAM=3000000 INPUT='$output' OUTPUT='$alldups' METRICS_FILE='$alldupsmetrics'");
system("mv $alldups $output");

sub runcmd
{
        my $cmd = shift @_;
        print STDERR "$cmd\n";
        system($cmd);
}
