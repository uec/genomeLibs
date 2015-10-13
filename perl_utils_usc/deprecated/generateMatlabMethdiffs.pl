#!/usr/bin/env perl

use File::Basename;
use File::Temp qw/ tempfile tempdir /;
use strict;
use Getopt::Long;

my $USAGE = "generateBareWigs.pl outputPrefix sample1REPLACECHROM.bam sample1REPLACECHROM.bam .. or ..\n" . 
"generateBareWigs.pl outputPrefix sample.bam .. or ..";
my $genome = "hg18";
my $resolution = 50;
my $low = 5.0;
my $high = 35.0;
GetOptions ('low=f',\$low, 'high=f',\$high, 'resolution=i',\$resolution, 'genome=s',\$genome) || die "$USAGE\n";
#print STDERR "doBare=${doBare}\tintervalFile=${intervalFile}\n";

# Input params
die "$USAGE\n" unless (@ARGV==3);
my ($newname , $csvTemplate1, $csvTemplate2) = @ARGV;

# Constants
my $MATLAB = "/usr/usc/matlab/2011a/bin/matlab";
my $TAG = "GTF2MATLABMETHDIFF";
my $GL = "/home/rcf-40/bberman/svn/genomeLibs/trunk";
my $classpath = "/home/rcf-40/bberman/Java/uecgatk/bin:/home/rcf-40/bberman/Java/uecgatk/lib/GenomeAnalysisTK.jar:/home/rcf-40/bberman/Java/uecgatk/lib/genomeLibs.jar:/home/rcf-40/bberman/Java/uecgatk/lib/StingUtils.jar:/home/rcf-40/bberman/Java/uecgatk/lib/UscKeck.jar:${GL}/sam-1.26.jar:${GL}/charts4j-1.2.jar:${GL}/biojava-live_1.6/apps-live.jar:${GL}/biojava-live_1.6/biojava-live.jar:${GL}/biojava-live_1.6/bytecode.jar:${GL}/commons-math-1.1.jar:${GL}/biojava-live_1.6/commons-cli.jar:${GL}/biojava-live_1.6/commons-collections-2.1.jar:${GL}/biojava-live_1.6/commons-dbcp-1.1.jar:${GL}/biojava-live_1.6/commons-pool-1.1.jar:${GL}/biojava-live_1.6/demos-live.jar:${GL}/biojava-live_1.6/jgrapht-jdk1.5.jar:${GL}/biojava-live_1.6/junit-4.4.jar:/home/uec-00/shared/production/software/picard/default/picard-tools-1.40.zip:${GL}/lib/tuple.jar";
my $matlabpath = "/home/uec-00/bberman/MatlabLibs/usc_matlab:/home/uec-00/bberman/MatlabLibs/insitu_matlab_g5:/home/uec-00/bberman/MatlabLibs/insitu_matlab";

my @FEATGTFS = ("~/tumor/genomic-data-misc/ENCODE/hg19/BroadHMM/BroadHMM-H1-hg19-StrongEnhancer.gtf");

# my $genomeToRef = {"hg18" => "hg18_unmasked", "hg19" => "hg19_rCRSchrm" };
# my $genomeToCpgsDir = {"hg18" => "Cpgls *.s", "hg19" => "CpgsHg19" };
# my $genomeRef = $genomeToRef->{$genome};
# die "Can't find genome ref for $genome\n" unless ($genomeRef);
# # Constants
# my $refDir = "/home/uec-00/shared/production/genomes/${genomeRef}";
# my $refFn = "${refDir}/${genomeRef}.fa";
my @CLUSTER_SIZES = qw/0 500/;

my @regions = ("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", 
	       "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", 
	       "chr19", "chr20", "chr21", "chr22", "chrX");
#@regions = ("chr2", "chr3", "chr5", "chr13");

@regions = ("chr4","chr11");
#my $CGINTERVALS_SEP = "/home/uec-02/bberman/BSseq/tumor/genomic-data-misc/Cpgs/${genomeRef}/CpgsRef.REPLACECHROM.7M-9M.bed";


#my $minCphFrac = 1.01;

my $lastChr = "";
my $lastEnd = 0;

# special case.  If newname contains a full path filename, remove it
$newname = basename($newname);

# Go through and check files
my $sepChrs = ($csvTemplate1 =~ /REPLACECHROM/);
foreach my $region (@regions)
{
    for my $csvTemplate ($csvTemplate1, $csvTemplate2)
    {
	my $fn = $csvTemplate;
	$fn =~ s/REPLACECHROM/${region}/g if ($sepChrs);
	die "Following file does not exist. Quitting\n${fn}\n" unless (-f $fn);
    }
}

# Start a fake command just to hold all jobs until the end
my $holdJobId = runCmd(0, "ls", "${TAG}_HOLDER", [], 1);
print STDERR "HOlding all jobs on job $holdJobId\n";


my @matlabJobs = ();
my $suffix = "";
#my @wigPrefixes = ();
#my @bedgraphPrefixes = ();
foreach my $region (@regions)
{
#     # 0-(nMinCsVariants) variable wig, #nMinCsVariants=bare wig
#     for (my $round = 0; $round <= $endround; $round++)
#     {
	my $csvFn1 = ${csvTemplate1};
	$csvFn1 =~ s/REPLACECHROM/${region}/g if ($sepChrs);
	my $base1 = basename($csvFn1);

	my $csvFn2 = ${csvTemplate2};
	$csvFn2 =~ s/REPLACECHROM/${region}/g if ($sepChrs);
	my $base2 = basename($csvFn2);
	my $prefix = "${newname}-${region}";
	my $cmd = "export CLASSPATH=${classpath}; export MATLABPATH=${matlabpath} ; ";

	my $lowToHighDesc = sprintf("low%.2f-high%.2f",$low,$high);
	$suffix = "-${lowToHighDesc}";
	$cmd .= "echo \"a=\'${csvFn1}\'; b=\'${csvFn2}\'; ";
	$cmd .= " processWGBSdiffMeth(a,b,\'$base1\',\'$base2\',\'${prefix}${suffix}\',$low,$high,[],0,$resolution);\" ";
	$cmd .= " | ${MATLAB} -nodisplay -nosplash -logfile ${prefix}-logger.txt";

	# It adds stuff on
	$suffix .= ".${resolution}bpResolution";
	my $gtfBase = "${newname}-${region}${suffix}";
	

	# Clustered
        # My stupid script leaves the last column blank, which doesn't work with IGV.  Also, i'm not sure why but chrX is showing
	# up as chr24 which becomes chrY.  Fix both
	my $fixGtfPipe = "perl -ne \'chomp; \@f=split(/\\t/);\$f[0]=~s\/chrY\/chrX\/g;\$f[8]=\"\"; print join(\"\\t\",\@f).\"\\n\";\'"; 
	foreach my $clusterSize (@CLUSTER_SIZES)
	{
	    $cmd .= "; java GtfClusteredWindows ${clusterSize} 0 0 0 none ${gtfBase}.gtf | $fixGtfPipe > ${gtfBase}.clustered${clusterSize}.gtf";
	}

	# And no need to keep original one
	$cmd .= "; rm -f ${gtfBase}.gtf ";
	
	my @dependJobs = ($holdJobId);
	my $curjobid = runCmd(0,$cmd, "${TAG}_matlab_r", \@dependJobs);
	push(@matlabJobs, $curjobid);
}

# Concatenate
my @catJobs = ();
my @catOuts = ();
if (scalar(@regions)>0)
{
    # No need to do concatenate "none" one
    # foreach my $clusterSize ("none", @CLUSTER_SIZES)
    foreach my $clusterSize (@CLUSTER_SIZES)
    {
	my $extension = "clustered${clusterSize}.gtf";
	$extension = "gtf" if ($clusterSize eq "none");
	
	my ($cmd, $fnout) = concatCmd(\@regions, $newname, $suffix, $extension);
	my @dependJobs = @matlabJobs;
	my $curjobid = runCmd(0,$cmd, "${TAG}_concat_clust${clusterSize}", \@dependJobs);
	push(@catJobs, $curjobid);
	push(@catOuts, $fnout);
    }
}

# Process each combined one
my @randOuts = ();
my @intersectionOuts = ();
foreach my $catOut (@catOuts)
{
    my @dependJobs = @catJobs;

    # And make a randomized control as well
    my $cmd .= "java edu.usc.epigenome.scripts.RandomizedIntervals -numTrials 1 -genome hg19 $catOut ";
    my $randOut = $catOut;
    $randOut =~ s/\.gtf$/.randomizedLocs.trial1.gtf/g;
    push(@randOuts, $randOut);

    foreach my $fnout ($catOut, $randOut)
    {
	my $intersectionOut = $fnout;
	$intersectionOut =~ s/\.gtf/.csv/g;
	my $featsSec = join(" ", @FEATGTFS); 
	$cmd .= "; java GtfFunctionalIntersection -genome hg19 $fnout $featsSec > $intersectionOut";
	push(@intersectionOuts, $intersectionOut);
    }

     my $curjobid = runCmd(0,$cmd, "${TAG}_postProcess", \@dependJobs);
}

# # Now igvtools and gzip
# foreach my $wigfn (@catOuts)
# {
#     my @dependJobs = @catJobs;
#     my $tdffn = $wigfn;
#     my @gzipDepends = @dependJobs;
#     if ($tdffn =~ /\.wig/)
#     {
# 	$tdffn =~ s/\.wig/\.tdf/g;
# 	my $cmd = "${IGVTOOLS} tile -c -f mean,min,max,median ${wigfn} ${tdffn} ${genome}";
# 	my $curjobid = runCmd(0,$cmd, "${TAG}_igvtools", \@dependJobs);
# 	push(@gzipDepends, $curjobid);
#     }

#     # Don't start the gzip until the tdf is finished, otherwise uncompressed file will disappear
#     my $cmd = "gzip ${wigfn}";
#     my $curjobid = runCmd(0,$cmd, "${TAG}_gzip", \@gzipDepends);
# }
# }

# Now we're ready to release the hounds
my $nsecs = 7;
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


# - - - - - - - - OLD 




# Returns jobid
sub runCmd
{
    my ($outputdir, $cmd, $prefix, $dependJobs, $userHold) = @_;

# I am sometimes getting "Could not create the Java virtual machine.", even with 7995mb.
    my ($fh, $file) = tempfile( "${prefix}XXXXXX" , DIR => "/tmp");
    print $fh "#Run on 1 processors on laird\n";
    print $fh "#PBS -l walltime=30:00:00\n#PBS -l mem=11995mb\n#PBS -l arch=x86_64\n#PBS -q lairdprio\n";
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

