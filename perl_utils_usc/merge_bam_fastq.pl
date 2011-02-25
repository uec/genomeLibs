#!/usr/bin/perl
use Getopt::Long;

#set program locations
$java = "/home/uec-00/shared/production/software/java/1.6.0_21/bin/java -Xmx2g";
$picard = "/home/uec-00/shared/production/software/picard/default";
$samtools = "/home/uec-00/shared/production/software/samtools/samtools";

#get cmd options
GetOptions(	'fastq=s' => \$inputFastq,
		'bam=s' => \$inputBam,
#		'readgroupid=s' => \$readgroupid,
		'samplename=s' => \$samplename,
		'libraryname=s' => \$libraryname,
#		'platform=s' => \$platform,
		'flowcell=s' => \$flowcell,
		'lane=s' => \$lane,
		'barcode=s' => \$barcode,
#		'platformunit=s' => \$platformunit,
#		'sequencingcenter=s' => \$sequencingcenter,
		'rundate=s' => \$rundate,
		'refgenome=s' => \$refGenome,
		'program=s' => \$program,
		'programversion=s' => \$programVersion,
		'programcmd=s' => \$programCmd,
		'ispaired=s' => \$isPaired,
		'isbisulfite=s' => \$isBisulfite,
		'output=s' => \$outputBam
				 	);
$readgroupid = $flowcell . "." . $barcode . $lane;

############################################################################
#convert fastq to bam
print STDERR "\nconvert reads to bam...\n";
$fastqbam = $inputFastq . ".bam";
runcmd("$java -jar $picard/FastqToSam.jar FASTQ='$inputFastq' QUALITY_FORMAT=Illumina OUTPUT='$fastqbam' READ_GROUP_NAME='$readgroupid' SAMPLE_NAME='$samplename' LIBRARY_NAME='$libraryname' PLATFORM_UNIT='$readgroupid' PLATFORM='illumina' SEQUENCING_CENTER='USC Epigenome Center' RUN_DATE='$rundate' SORT_ORDER='queryname'");

############################################################################
#prepare aligned reads:
print STDERR "\npreprocess the aligned reads...\n";
$alnsam = $inputBam . ".sam";
$alnsamclean = $inputBam . ".clean.sam";

unlink $alnsam;
unlink $alnsamclean;

runcmd("$samtools view -h $inputBam -o $alnsam");

#fix maq mangling
$queryname = `head -n 1 $inputFastq`;
$queryname =~ /\@(.+?)\:/;
$queryname = $1;
print STDERR "fixing maq-mangled querynames: $queryname\n";


open(ALNSAM, "<$alnsam");
open(ALNSAMCLEAN, ">$alnsamclean");
print ALNSAMCLEAN `$samtools view -H $fastqbam`;
while($line = <ALNSAM>)
{
	$line =~ s/^.+?\:/$queryname\:/ if($line !~ /^\@/);
	print ALNSAMCLEAN $line;	
}
close ALNSAM;
close ALNSAMCLEAN;

#add readgroups
print STDERR "fixing read groups in aligned bam\n";
runcmd("sed -i '/*/s|\$|\\tRG:Z:'$readgroupid'|' $alnsamclean");

#######################################################################
#merge aln and unaln
print STDERR "\nrecombine aln and un-aln...\n";
runcmd("$java -jar $picard/MergeBamAlignment.jar UNMAPPED_BAM='$fastqbam' ALIGNED_BAM='$alnsamclean' OUTPUT='$outputBam' REFERENCE_SEQUENCE='$refGenome' PROGRAM_RECORD_ID='$program' PROGRAM_GROUP_VERSION='$programVersion' PROGRAM_GROUP_COMMAND_LINE='$programCmd' PROGRAM_GROUP_NAME='$program' PAIRED_RUN=$isPaired IS_BISULFITE_SEQUENCE=$isBisulfite ALIGNED_READS_ONLY=false");

#######################################################################
#mark duplicates
print STDERR "\nmark dups...\n";
$alldups = $outputBam . ".dups.bam";
$alldupsmetrics = $outputBam . ".dups.bam.metrics";
runcmd("$java -jar $picard/MarkDuplicates.jar INPUT='$outputBam' OUTPUT='$alldups' METRICS_FILE='$alldupsmetrics'");

#######################################################################
#clean up
print STDERR "\ncleaning up...\n";
runcmd("mv $alldups $outputBam");
unlink $alnsam;
unlink $alnsamclean;
unlink $fastqbam;

exit;

sub runcmd
{
	my $cmd = shift @_;
	print STDERR "$cmd\n";
	system($cmd);
}
