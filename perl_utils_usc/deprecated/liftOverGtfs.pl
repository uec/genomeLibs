#!/usr/bin/perl

use strict;
use File::Basename;
use File::Temp qw/ tempfile tempdir /;
use strict;
use Getopt::Long;

my $TAG = "LIFTOVER";
my $EXEC = "/home/uec-00/bberman/bin/liftOver";
my $CHAINDIR = "/home/uec-00/bberman/lib";
my $MEMPERJOB = 22995; # 7000  # The wig_to_bedgraph is ridiculously inefficient takes up to 16G or more. Should rewrite

my $USAGE = "liftOverGtfs.pl hg18(from) hg19(to) file1.wig file2.wig ... (if hg18=hg18, no liftover)";

die "$USAGE\n" if (@ARGV < 3);

# GetOptions ('MEMPERJOB=i' => \$MEMPERJOB, 'doGnome!' => \$doGnome, 'minMapQ=i' => \$minMapQ,  
# 'minCT=i' => \$minCT, 'cpus=i' => \$cpus, 'minContextFracReadsMatching=f' => \$minContextFracReadsMatching,
# 'minNumWcg=i' => \$minNumWcg, 'minNumGch=i' => \$minNumGch, 'contextCombos!' => \$contextCombos,
# 'readBased!' => \$readBased, 'hcgs!' => \$doHcgs, 'ccgs!' => \$doCcgs) || die "$USAGE\n";


my $fromVers = shift(@ARGV);
my $toVers = shift(@ARGV);
my (@fns) = @ARGV;

# Start a fake command just to hold all jobs until the end
my $holdJobId = runCmd(0, "ls", "${TAG}_HOLDER", [], 1);
print STDERR "HOlding all jobs on job $holdJobId\n";

my @jobs = ();
GTF: foreach my $fn (@fns)
{
    my $cmd = "ls"; # just so our next semicolon is valid

    # Need a temp file for unmapped
    my $tmpfn = "/tmp/liftoverrand." . int(rand(100000000));

    # Make a base file name
    my $basefn = $fn;
    $basefn =~ s/\.g[tf]f$//gi;
    $basefn =~ s/\.?${fromVers}$//gi;

    my $fnToVers;
    if ($fromVers eq $toVers)
    {
	$fnToVers = $fn;
    }
    else
    {
	# First , get rid of trackline
	my $chrSec = "";
#	$chrSec = " | grep chr3 ";  ##### ------ ##### ----- TESTING ----- ####### --- 
	$cmd .= " ; grep -v 'track' ${fn} ${chrSec} > ${fn}.temp ";

	# Now do liftover of bedgraph
	$fnToVers = "${basefn}.${toVers}.gtf";

	# Get chain file
	my $toVersChain = $toVers;
	substr($toVersChain,0,1,uc(substr($toVersChain,0,1)));
	my $chainBare = "${fromVers}To${toVersChain}.over.chain";
	my $chainfn = "${CHAINDIR}/${chainBare}";
	if (!(-f $chainfn))
	{
	    print STDERR "Can't read chain file \"$chainfn\"\n";
	    last GTF; # None of them will work.  But we don't want to die because we need to release holder job.
	}
	$cmd .= "; ${EXEC} -gff ${fn}.temp ${chainfn} ${fnToVers} ${tmpfn}";
    }

    # Cleanup
    $cmd .= " ; rm -f $tmpfn ${fn}.temp";

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
