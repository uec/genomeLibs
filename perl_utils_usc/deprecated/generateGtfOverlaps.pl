#!/usr/bin/perl

use strict;
use File::Basename;
use File::Temp qw/ tempfile tempdir /;
use strict;
use Getopt::Long;

my $TAG = "GTFOVS";
my $EXEC = "/home/uec-00/bberman/bin/liftOver";
my $MEMPERJOB = 7500;
my $MINMEM_JAVA = $MEMPERJOB - 2500; # Above this and Java can fail to be launched

my $USAGE = "generateGtfOverlaps.pl promoters.gtf promoters-pretty-name intersect1.gtf intersect2.gtf ...";

die "$USAGE\n" if (@ARGV < 3);

# GetOptions ('MEMPERJOB=i' => \$MEMPERJOB, 'doGnome!' => \$doGnome, 'minMapQ=i' => \$minMapQ,  
# 'minCT=i' => \$minCT, 'cpus=i' => \$cpus, 'minContextFracReadsMatching=f' => \$minContextFracReadsMatching,
# 'minNumWcg=i' => \$minNumWcg, 'minNumGch=i' => \$minNumGch, 'contextCombos!' => \$contextCombos,
# 'readBased!' => \$readBased, 'hcgs!' => \$doHcgs, 'ccgs!' => \$doCcgs) || die "$USAGE\n";

my $classpath = "/home/rcf-40/bberman/Java/uecgatk/bin:/home/rcf-40/bberman/Java/uecgatk/lib/GenomeAnalysisTK.jar:/home/rcf-40/bberman/Java/uecgatk/lib/tuple.jar:/home/rcf-40/bberman/Java/uecgatk/bin:/home/rcf-40/bberman/Java/uecgatk/lib/GenomeAnalysisTK.jar:/home/rcf-40/bberman/Java/uecgatk/lib/genomeLibs.jar:/home/rcf-40/bberman/Java/uecgatk/lib/StingUtils.jar:/home/rcf-40/bberman/Java/uecgatk/lib/UscKeck.jar:/home/rcf-40/bberman/svn/genomeLibs/trunk/sam-1.26.jar:/home/rcf-40/bberman/svn/genomeLibs/trunk/charts4j-1.2.jar:/home/rcf-40/bberman/svn/genomeLibs/trunk/biojava-live_1.6/apps-live.jar:/home/rcf-40/bberman/svn/genomeLibs/trunk/biojava-live_1.6/biojava-live.jar:/home/rcf-40/bberman/svn/genomeLibs/trunk/biojava-live_1.6/bytecode.jar:/home/rcf-40/bberman/svn/genomeLibs/trunk/commons-math-1.1.jar:/home/rcf-40/bberman/svn/genomeLibs/trunk/biojava-live_1.6/commons-cli.jar:/home/rcf-40/bberman/svn/genomeLibs/trunk/biojava-live_1.6/commons-collections-2.1.jar:/home/rcf-40/bberman/svn/genomeLibs/trunk/biojava-live_1.6/commons-dbcp-1.1.jar:/home/rcf-40/bberman/svn/genomeLibs/trunk/biojava-live_1.6/commons-pool-1.1.jar:/home/rcf-40/bberman/svn/genomeLibs/trunk/biojava-live_1.6/demos-live.jar:/home/rcf-40/bberman/svn/genomeLibs/trunk/biojava-live_1.6/jgrapht-jdk1.5.jar:/home/rcf-40/bberman/svn/genomeLibs/trunk/biojava-live_1.6/junit-4.4.jar:/home/uec-00/shared/production/software/picard/default/picard-tools-1.40.zip";

my ($targfn, $targname, @intfns) = @ARGV;

# Make sure the targfn is a file and targname isn't
die "Can't find file $targfn\n$USAGE\n" if (!(-e $targfn));
die "Why is 2nd argument a file?\n$USAGE\n" if (-e $targname);

# Start a fake command just to hold all jobs until the end
my $holdJobId = runCmd(0, "ls", "${TAG}_HOLDER", [], 1);
print STDERR "HOlding all jobs on job $holdJobId\n";

my @jobs = ();
GTF: foreach my $intfn (@intfns)
{
    my $cmd = "ls"; # just so our next semicolon is valid

    # If it's gzipped, unzip
    my $gz = 0;
    if ($intfn =~ /\.gz$/)
    {
	$gz = 1;
	$cmd .= " ; gunzip $intfn";
	$intfn =~ s/\.gz$//g;
    }

    # If it's a peaks file, we have to translate it
    my $gtffn = $intfn;
    if ($intfn =~ /\.peaks$/)
    {
	$cmd .= " ; encodePeakFilesToGtf.pl $intfn";
	my ($base, $path, $suf) = fileparse($intfn, qr/\.[^.]*/);
	$gtffn = $base . ".peaks.gtf";
    }

    # call into java
    my ($base, $path, $suf) = fileparse($gtffn, qr/.g[tf]f/i);
    my $outname = "${targname}-${base}.gtf";
    my $exportSec = "export CLASSPATH=${classpath}";
    $cmd .= "; $exportSec ; java -Xmx${MINMEM_JAVA}m  GtfFilterOverlapping 0 ${targfn} ${gtffn} > $outname ";

    # Cleanup
    $cmd .= " ; rm -f $gtffn"  if ($gtffn ne $intfn);

    my @dependJobs = ($holdJobId);
    
    my $curjobid = runCmd(0,$cmd, "${TAG}_main", \@dependJobs);
    push(@jobs, $curjobid);
}


# Now we're ready to release the hounds
my $nsecs = 1;
print STDERR "About to start jobs by releasing $holdJobId after a ${nsecs} second break\n";
sleep($nsecs);
`qrls $holdJobId`;



# - - - - - - - - - 
# Returns jobid
sub runCmd
{
    my ($outputdir, $cmd, $prefix, $dependJobs, $userHold) = @_;

# I am sometimes getting "Could not create the Java virtual machine.", even with 7995mb.
    my ($fh, $file) = tempfile( "${prefix}XXXXXX" , DIR => "/tmp");
    print $fh "#Run on 1 processors on laird\n";
    my $minmem = $MEMPERJOB+500;
    # When HPCC 10gb is down
    print $fh "#PBS -l walltime=30:00:00\n#PBS -l mem=${minmem}mb\n#PBS -l arch=x86_64\n#PBS -q lairdprio\n";
#    print $fh "#PBS -l walltime=00:30:00\n#PBS -l mem=${minmem}mb\n";
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
