#!/usr/bin/env perl

use File::Basename;
use File::Temp qw/ tempfile tempdir /;
use strict;
use Getopt::Long;

my $TAG = "BBC2REG";

my $USAGE = "generateBamsByChromToMethylChromTables.pl [--noremoveTemps] [--useReadCGschema] [--dumpCoverageMatrices] tablePrefix sampleREPLACECHROM.bam samplexREPLACEHROM.bam ...";
my @regions = ("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", 
	       "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", 
	       "chr19", "chr20", "chr21", "chr22", "chrX", "chrY", "chrM");
#@regions = ("chr20","chr21");
#@regions = ("chrM");
#@regions = ("chr2", "chr3", "chr5", "chr13");

my $SAMDIR = "/home/uec-00/shared/production/software/samtools";
my $RMTMPS = 0;
my $useReadCGschema = 0;
my $minCphFrac = 1.01;
my $minCphCoverage = 0;
my $dumpCoverageMatrices = 0;
GetOptions ('removeTemps!' => \$RMTMPS, 'useReadCGschema!' => \$useReadCGschema, 'dumpCoverageMatrices!' => \$dumpCoverageMatrices, 'minCphFrac' => \$minCphFrac) || die "$USAGE\n";

print STDERR "removeTemps=${RMTMPS}\n";

# Input params
die "$USAGE\n" unless (@ARGV>=2);
my ($newname , @bamTemplates) = @ARGV;

# Go through and check files
foreach my $bamTemplate (@bamTemplates)
{
    die "Following file does not have the string \"REPLACECHROM\" . Quitting\n${bamTemplate}\n" unless ($bamTemplate =~ /REPLACECHROM/);
    foreach my $region (@regions)
    {
	my $fn = $bamTemplate;
	$fn =~ s/REPLACECHROM/${region}/g;
	die "Following file does not exist. Quitting\n${fn}\n" unless (-f $fn);
    }
}


# Start a fake command just to hold all jobs until the end
my $holdJobId = runCmd(0, "ls", "${TAG}_HOLDER", [], 1);
print STDERR "HOlding all jobs on job $holdJobId\n";

foreach my $region (@regions)
{
    my $tmpDir = "./tmp_${newname}_${region}";
    `mkdir $tmpDir`;
    my @dependJobs = ($holdJobId);
    my @mergedBamFns = ();

    # Merge maps within region directory
    my $mergeBamFn = "";
    my $mergeJobids = [];
    # If there's only 1, no need to merge
    if (scalar(@bamTemplates)<2)
    {
	$mergeBamFn = $bamTemplates[0];
	$mergeBamFn =~ s/REPLACECHROM/${region}/g;
	push(@$mergeJobids, @dependJobs);
    }
    else
    {
	my @inputBams = ();
	foreach my $bamTemplate (@bamTemplates)
	{
	    my $fn = $bamTemplate;
	    $fn =~ s/REPLACECHROM/${region}/g;
	    push(@inputBams, $fn);
	}
	
	my $mergeName = "${tmpDir}/${newname}_${region}";
	my $mergeJobid;
	($mergeJobid, $mergeBamFn) = mergeBams($mergeName, \@dependJobs, \@inputBams, $tmpDir);
	push(@$mergeJobids, $mergeJobid);
    }

    # Convert to methyldb table
    my $cmd = "export CLASSPATH=/home/rcf-40/bberman/svn/genomeLibs/trunk/genomeLibs.jar:/home/rcf-40/bberman/svn/genomeLibs/trunk/UscKeck.jar:/home/rcf-40/bberman/svn/genomeLibs/trunk/sam-1.26.jar:/home/rcf-40/bberman/svn/genomeLibs/trunk/biojava-live_1.6/apps-live.jar:/home/rcf-40/bberman/svn/genomeLibs/trunk/biojava-live_1.6/biojava-live.jar:/home/rcf-40/bberman/svn/genomeLibs/trunk/biojava-live_1.6/bytecode.jar:/home/rcf-40/bberman/svn/genomeLibs/trunk/commons-math-1.1.jar:/home/rcf-40/bberman/svn/genomeLibs/trunk/biojava-live_1.6/commons-cli.jar:/home/rcf-40/bberman/svn/genomeLibs/trunk/biojava-live_1.6/commons-collections-2.1.jar:/home/rcf-40/bberman/svn/genomeLibs/trunk/biojava-live_1.6/commons-dbcp-1.1.jar:/home/rcf-40/bberman/svn/genomeLibs/trunk/biojava-live_1.6/commons-pool-1.1.jar:/home/rcf-40/bberman/svn/genomeLibs/trunk/biojava-live_1.6/demos-live.jar:/home/rcf-40/bberman/svn/genomeLibs/trunk/biojava-live_1.6/jgrapht-jdk1.5.jar:/home/rcf-40/bberman/svn/genomeLibs/trunk/biojava-live_1.6/junit-4.4.jar\n";

    if ($useReadCGschema)
    {
    	$cmd .= "java -Xmx3500m edu.usc.epigenome.scripts.SamToMethylreadCGdbOffline -chrom ${region} -minMapQ 30 ${newname} $mergeBamFn";
    }
    elsif ($dumpCoverageMatrices)
    {
	$cmd .= "java -Xmx3500m edu.usc.epigenome.scripts.SamToConversionByCoverageMatrix -chrom ${region} -minMapQ 30 -outputCphs -minNextBaseCoverage 5 -minOppStrandCoverage 5 -maxCoverageOutput 300 ${newname} $mergeBamFn";
    }
    else
    {
	$cmd .= "java -Xmx3500m edu.usc.epigenome.scripts.SamToMethyldbOffline -minConv 1 -useCpgsToFilter -minMapQ 30 -chrom ${region} ${newname} $mergeBamFn";
    }

    my $curJobids = [runCmd($tmpDir,$cmd, "${TAG}_createMehylTable", $mergeJobids)];

}

# Now we're ready to release the hounds
my $nsecs = 20;
print STDERR "About to start jobs by releasing $holdJobId after a ${nsecs} second break\n";
sleep($nsecs);
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
	$curJobids = [runCmd($tmpdir, $cmd, "${TAG}_mvBams", $dependJobs)];
    }
    else
    {
	my $cmd = "${SAMDIR}/samtools merge ${curOut} " . join(" ", @$individualBamFns);
	$cmd .= "; rm -f " . join(" ", @$individualBamFns) if ($RMTMPS);
	$curJobids = [runCmd($tmpdir, $cmd, "${TAG}_mergeBams", $dependJobs)];
    }

    # Do not remove dups
    my $finalBamOut = $curOut;

    # Index new bam
    my $curIn = $curOut;
    $curOut = "${curIn}.bai";
    my $cmd = "${SAMDIR}/samtools index ${curIn} ${curOut}";
    $curJobids = [runCmd($tmpdir, $cmd, "${TAG}_index", $curJobids)];

    return (@${curJobids}[0], $finalBamOut);
}




# - - - - - - - - OLD 




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
	`echo "job ${lineout} =" >> ${TAG}_pbs_log.txt`;
    }

    `cat $file >> ${TAG}_pbs_log.txt`;

    unlink($file);
    return $lineout;
}

