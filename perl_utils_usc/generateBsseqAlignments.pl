#!/usr/bin/env perl

use File::Basename;
use File::Temp qw/ tempfile tempdir /;
use strict;
use Getopt::Long;

my $TAG = "BSSEQ_ALIGN";
#my @SUFFIXES = qw/GCH HCG readcvg/;


my $USAGE = "generateBsseqAlignments.pl flankSize downsampleCols bsseq1.bam bsseq2.bam , file1.gff file2.gff ...";

my $classpath = "/home/rcf-40/bberman/Java/uecgatk/bin:/home/rcf-40/bberman/Java/uecgatk/lib/GenomeAnalysisTK.jar:/home/rcf-40/bberman/Java/uecgatk/lib/tuple.jar:/home/rcf-40/bberman/Java/uecgatk/bin:/home/rcf-40/bberman/Java/uecgatk/lib/GenomeAnalysisTK.jar:/home/rcf-40/bberman/Java/uecgatk/lib/genomeLibs.jar:/home/rcf-40/bberman/Java/uecgatk/lib/StingUtils.jar:/home/rcf-40/bberman/Java/uecgatk/lib/UscKeck.jar:/home/rcf-40/bberman/svn/genomeLibs/trunk/sam-1.26.jar:/home/rcf-40/bberman/svn/genomeLibs/trunk/charts4j-1.2.jar:/home/rcf-40/bberman/svn/genomeLibs/trunk/biojava-live_1.6/apps-live.jar:/home/rcf-40/bberman/svn/genomeLibs/trunk/biojava-live_1.6/biojava-live.jar:/home/rcf-40/bberman/svn/genomeLibs/trunk/biojava-live_1.6/bytecode.jar:/home/rcf-40/bberman/svn/genomeLibs/trunk/commons-math-1.1.jar:/home/rcf-40/bberman/svn/genomeLibs/trunk/biojava-live_1.6/commons-cli.jar:/home/rcf-40/bberman/svn/genomeLibs/trunk/biojava-live_1.6/commons-collections-2.1.jar:/home/rcf-40/bberman/svn/genomeLibs/trunk/biojava-live_1.6/commons-dbcp-1.1.jar:/home/rcf-40/bberman/svn/genomeLibs/trunk/biojava-live_1.6/commons-pool-1.1.jar:/home/rcf-40/bberman/svn/genomeLibs/trunk/biojava-live_1.6/demos-live.jar:/home/rcf-40/bberman/svn/genomeLibs/trunk/biojava-live_1.6/jgrapht-jdk1.5.jar:/home/rcf-40/bberman/svn/genomeLibs/trunk/biojava-live_1.6/junit-4.4.jar:/home/uec-00/shared/production/software/picard/default/picard-tools-1.40.zip";

my $refFn = "/home/uec-00/bberman/genomes/hg18_unmasked/hg18_unmasked.plusContam.fa";
my $gffToBedBin = "/home/uec-00/bberman/svn/genomeLibs/trunk/perl_utils_usc/gffToBed.pl";

#my $minCphFrac = 1.01;
my $minCT = 3;
my $doGnome = 1;
my $minMapQ = 30;
my $MEMPERJOB = 2000;
GetOptions ('MEMPERJOB=i', \$MEMPERJOB, 'doGnome!' => \$doGnome, 'minMapQ=i',\$minMapQ,  'minCT=i',\$minMapQ) || die "$USAGE\n";

# Input params
die "$USAGE\n" unless (@ARGV>=5);
my $flank = shift(@ARGV);
my $downscaleCols = shift(@ARGV);
my $commaSeen = 0;
my $bams = [];
my $gffs = [];
foreach my $arg (@ARGV)
{
    if ($arg eq ',')
    {
	$commaSeen = 1;
	next;
    }

    my $l = ($commaSeen) ? $gffs : $bams;
    push(@$l, $arg);
}
print STDERR sprintf("Found %d bams and %d gffs\n",scalar(@$bams), scalar(@$gffs));

# Start a fake command just to hold all jobs until the end
my $holdJobId = runCmd(0, "ls", "${TAG}_HOLDER", [], 1);
print STDERR "HOlding all jobs on job $holdJobId\n";

my @jobs = ();

foreach my $bam (@$bams)
{
    foreach my $gff (@$gffs)
    {
	my $exportSec = "export CLASSPATH=${classpath}";

	# Make an interval file, which has to be BED
#	my $cmd = " $exportSec ; java -Xmx${MEMPERJOB}m  GtfFlankedWindows ${flank} ${gff}"; # This is not replacing chr24 correctly
	my $flankPlus = int($flank*1.3);
	my $cmd = " $exportSec ; java -Xmx${MEMPERJOB}m  GtfCenteredWindows ${flankPlus} ${gff}";

	# Has to be bed
	my ($bamBase) = basename($bam);
	my ($gffBase) = basename($gff);
	my $flankedGff = sprintf("%s.%dbp.gtf",$gffBase,$flank);
	my $flankedBed = sprintf("%s.%dbp.bed",$gffBase,$flank);
	$cmd .= " ; ${gffToBedBin} $flankedGff";

	my $cphSec = ($doGnome) ? " -cph " : "";
	my $outfn = sprintf("ALIGN_%d_%d_mapq%d_minct%d.%s.%s.csv", $flank, $downscaleCols, $minMapQ, $minCT, $bamBase, $gffBase);
	$cmd .= " ; java -Xmx${MEMPERJOB}m org.broadinstitute.sting.gatk.CommandLineGATK -T  GnomeSeqFeatureAlignment -R ${refFn} -I ${bam} --min_mapping_quality_score ${minMapQ} ${cphSec} --minCT ${minCT} --intervals ${flankedBed} --elementGff ${gff} --windSize ${flank} -et NO_ET -nt 1 --downscaleCols ${downscaleCols} --outPrefix ${outfn} ";

# 	# Clean up
# 	$cmd .= " ; rm -f ${flankedGff} ${flankedBed}";

	my @dependJobs = ($holdJobId);
    
	my $curjobid = runCmd(0,$cmd, "${TAG}_main", \@dependJobs);
	push(@jobs, $curjobid);
    }
}

# Now we're ready to release the hounds
my $nsecs = 3;
print STDERR "About to start jobs by releasing $holdJobId after a ${nsecs} second break\n";
sleep($nsecs);
`qrls $holdJobId`;



# - - - - - - - -

sub concatCmd
{
    my ($regions, $prefix, $suffix, $extension) = @_;

    my $fnout = "${prefix}${suffix}.${extension}";
    my $tmpfnout = "tmp${fnout}";

    my $cmd = "";
    my $count = 0;
    foreach my $region (@$regions)
    {
	my $fn = "${prefix}-${region}${suffix}.${extension}";
	$cmd .= ";\n echo \"\" >> $tmpfnout" if ($count>0); # They don't put a newline
	my $grepsec = ($count == 0) ? "" : " | grep -v track ";
	my $arrow = ($count == 0) ? ">" : ">>";
	$cmd .= ";\n " if ($count>0);
	$cmd .= "cat $fn ${grepsec} ${arrow} ${tmpfnout}";

	$count++;
    }

    # The newlines are unpredictable. Get rid of bare lines
    $cmd .= ";\n egrep -v \'\^\\s\*\$\' ${tmpfnout} > $fnout;\n rm -f $tmpfnout";

    return ($cmd, $fnout);
}


# Returns jobid
sub runCmd
{
    my ($outputdir, $cmd, $prefix, $dependJobs, $userHold) = @_;

# I am sometimes getting "Could not create the Java virtual machine.", even with 7995mb.
    my ($fh, $file) = tempfile( "${prefix}XXXXXX" , DIR => "/tmp");
    print $fh "#Run on 1 processors on laird\n";
    my $minmem = $MEMPERJOB+500;
    print $fh "#PBS -l walltime=30:00:00\n#PBS -l mem=${minmem}mb\n#PBS -l arch=x86_64\n#PBS -q lairdprio\n";
    # PBS can crash with too many afterany on the same line
    print $fh "#PBS -W depend=afterany:" . join(":",@$dependJobs) ."\n" if ($dependJobs && @$dependJobs);
#     if ($dependJobs)
#     {
# 	foreach my $dj (@$dependJobs)
# 	{
# 	    print $fh "#PBS -W depend=afterany:$dj\n";
# 	}
#     }
    print $fh "cd \"\$PBS_O_WORKDIR\"\n";
    print $fh "${cmd}\n\n";
    close($fh);

    my $outputSec = ($outputdir) ? " -e ${outputdir} -o ${outputdir} " : "";
    my $holdSec = ($userHold) ? " -h " : "";
    my $fullcmd = "qsub ${outputSec} ${holdSec} $file";
    print STDERR "${cmd}\n${fullcmd}\n";

    my $lineout = "";
    if (1)
    {
	$lineout = `$fullcmd`;
	chomp $lineout;
	print STDERR "runCmd got lineout:\t\"$lineout\"\n";
	`echo "job ${lineout} =" >> ${TAG}_pbs_log.txt`;
    }

    `cat $file >> ${TAG}_pbs_log.txt`;

    unlink($file);
    return $lineout;
}

