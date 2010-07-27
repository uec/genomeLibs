#!/usr/bin/env perl

use File::Basename;
use File::Temp qw/ tempfile tempdir /;
use strict;
use Getopt::Long;

my $USAGE = "generateBamsToMethylChromTables.pl [--noremoveTemps] [--useReadCGschema] tablePrefix f1.bam f2.bam f3.bam ...";
my @regions = ("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", 
	       "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", 
	       "chr19", "chr20", "chr21", "chr22", "chrX", "chrY", "chrM");
#@regions = ("chr20","chr21");

my $SAMDIR = "/home/uec-00/bberman/bin";
my $RMTMPS = 0;
my $useReadCGschema = 0;
GetOptions ('removeTemps!' => \$RMTMPS, 'useReadCGschema!' => \$useReadCGschema) || die "$USAGE\n";

print STDERR "removeTemps=${RMTMPS}\n";

# Input params
die "$USAGE\n" unless (@ARGV>=2);
my ($newname , @bamFns) = @ARGV;

# Check if we have any postprocessing steps to add to filenames
my $firstbam = $bamFns[0];
my $processingSec = "";
if ($firstbam =~ /(\.NODUPS.*)\.bam$/)
{
    $processingSec = $1;
}
print STDERR "proc=${processingSec}\n";

# Start a fake command just to hold all jobs until the end
my $holdJobId = runCmd(0, "ls", "B2REG_HOLDER", [], 1);
print STDERR "HOlding all jobs on job $holdJobId\n";

foreach my $region (@regions)
{
    my $tmpDir = "./tmp_${newname}_${region}";
    `mkdir $tmpDir`;
    my @dependJobs = ($holdJobId);
    my @mergedBamFns = ();

    # Generate regional bams within region directory
    my @regionalBamFns = ();
    foreach my $bamFn (@bamFns)
    {
	if ($bamFn !~ /\.bam$/)
	{
	    print STDERR "File $bamFn does not end in .bam\n";
	}

	my ($name, $path, $suf) = fileparse($bamFn, qr/\.bam/);
	my $regionalBamFn = "${tmpDir}/${name}_${region}.bam";
	my $pullRegionJobid = pullRegions($bamFn, $regionalBamFn, $region, $tmpDir,\@dependJobs);
	push(@dependJobs, $pullRegionJobid);
	push(@regionalBamFns, $regionalBamFn);
    }

    # Merge maps within region directory
    my $mergeName = "${tmpDir}/${newname}_${region}${processingSec}";
    my ($mergeJobid, $mergeBamFn) = mergeBams($mergeName, \@dependJobs, \@regionalBamFns, $tmpDir);

    # Convert to methyldb table
    my $cmd = "export CLASSPATH=/home/rcf-40/bberman/svn/genomeLibs/trunk/genomeLibs.jar:/home/rcf-40/bberman/svn/genomeLibs/trunk/UscKeck.jar:/home/rcf-40/bberman/svn/genomeLibs/trunk/sam-1.07.jar:/home/rcf-40/bberman/svn/genomeLibs/trunk/biojava-live_1.6/apps-live.jar:/home/rcf-40/bberman/svn/genomeLibs/trunk/biojava-live_1.6/biojava-live.jar:/home/rcf-40/bberman/svn/genomeLibs/trunk/biojava-live_1.6/bytecode.jar:/home/rcf-40/bberman/svn/genomeLibs/trunk/commons-math-1.1.jar:/home/rcf-40/bberman/svn/genomeLibs/trunk/biojava-live_1.6/commons-cli.jar:/home/rcf-40/bberman/svn/genomeLibs/trunk/biojava-live_1.6/commons-collections-2.1.jar:/home/rcf-40/bberman/svn/genomeLibs/trunk/biojava-live_1.6/commons-dbcp-1.1.jar:/home/rcf-40/bberman/svn/genomeLibs/trunk/biojava-live_1.6/commons-pool-1.1.jar:/home/rcf-40/bberman/svn/genomeLibs/trunk/biojava-live_1.6/demos-live.jar:/home/rcf-40/bberman/svn/genomeLibs/trunk/biojava-live_1.6/jgrapht-jdk1.5.jar:/home/rcf-40/bberman/svn/genomeLibs/trunk/biojava-live_1.6/junit-4.4.jar\n";

	if ($useReadCGschema)
	{
    	$cmd .= "java -Xmx4000m edu.usc.epigenome.scripts.SamToMethylreadCGdbOffline -chrom ${region} -minMapQ 30 ${newname} $mergeBamFn";
	}
	else
	{
    	$cmd .= "java -Xmx4000m edu.usc.epigenome.scripts.SamToMethyldbOffline -minConv 1 -useCpgsToFilter -minMapQ 30 -chrom ${region} ${newname} $mergeBamFn";
	}

    my $curJobids = [runCmd($tmpDir,$cmd, "B2REG_createMehylTable", [$mergeJobid])];

}

# Now we're ready to release the hounds
print STDERR "About to start jobs by releasing $holdJobId after a 60 second break\n";
sleep(60);
`qrls $holdJobId`;



# Returns
# ($processId, $outBamFn)
sub mergeBams
{
    my ($prefix, $dependJobs, $individualBamFns, $tmpdir) = @_;

    print STDERR "merging bams (".join(",",@$individualBamFns) .")\n\tFrom jobs(".join(",",@$dependJobs).")\n";

    my $curOut = "${prefix}.bam";
    my $curJobids;
    if (scalar(@$individualBamFns) == 1)
    {
	my $cmd = "mv " . @{$individualBamFns}[0] . " $curOut";
	$curJobids = [runCmd($tmpdir, $cmd, "B2REG_mvBams", $dependJobs)];
    }
    else
    {
	my $cmd = "${SAMDIR}/samtools merge ${curOut} " . join(" ", @$individualBamFns);
	$cmd .= "; rm -f " . join(" ", @$individualBamFns) if ($RMTMPS);
	$curJobids = [runCmd($tmpdir, $cmd, "B2REG_mergeBams", $dependJobs)];
    }

    # We have to remove dups again
    my $curIn = $curOut;
    $curOut = "${prefix}.NODUPS.bam";
    my $cmd = "${SAMDIR}/samtools rmdup -s ${curIn} ${curOut}";
    $cmd .= "; rm -f ${curIn}" if ($RMTMPS);
    $curJobids = [runCmd($tmpdir,$cmd, "B2REG_rmdups", $curJobids)];

    my $finalBamOut = $curOut;

    # Index bam
    my $curIn = $curOut;
    $curOut = "${curIn}.bai";
    my $cmd = "${SAMDIR}/samtools index ${curIn} ${curOut}";
    $curJobids = [runCmd($tmpdir, $cmd, "B2REG_index", $curJobids)];

    return (@${curJobids}[0], $finalBamOut);
}




# - - - - - - - - OLD 




# Returns the final process ID
sub pullRegions
{
    my ($inBamFn, $outBamFn, $region, $tmpdir, $holdJobsIds) = @_;

    my $finalProcId = 0;

    if ($inBamFn !~ /\.bam$/)
    {
	print STDERR "pipeline only works with .bam input files!\n";
	exit;
    }

    my $curIn = $inBamFn;
    my $curOut = $outBamFn;
    my $curJobids = $holdJobsIds;

    my $cmd = "${SAMDIR}/samtools view -b -o ${curOut} ${curIn} ${region}";
    $curJobids = [runCmd($tmpdir,$cmd, "B2REG_pullRegion", $curJobids)];

    # Index bam
    $curIn = $curOut;
    $curOut = "${curOut}.bai";
    my $cmd = "${SAMDIR}/samtools index ${curIn} ${curOut}";
    $curJobids = [runCmd($tmpdir,$cmd, "B2REG_pullRegionIndex", $curJobids)];

    $finalProcId = @${curJobids}[0];
    return $finalProcId;
}

# Returns jobid
sub runCmd
{
    my ($outputdir, $cmd, $prefix, $dependJobs, $userHold) = @_;

    my ($fh, $file) = tempfile( "${prefix}XXXXXX" , DIR => "/tmp");
    print $fh "#Run on 1 processors on laird\n";
    print $fh "#PBS -l walltime=30:00:00\n#PBS -l mem=3995mb\n#PBS -l arch=x86_64\n#PBS -q laird\n";
    print $fh "#PBS -W depend=afterany:" . join(":",@$dependJobs) ."\n" if ($dependJobs && @$dependJobs);
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
	`echo "job ${lineout} =" >> B2REG_pbs_log.txt`;
    }

    `cat $file >> B2REG_pbs_log.txt`;

    unlink($file);
    return $lineout;
}

