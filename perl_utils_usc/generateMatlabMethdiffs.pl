#!/usr/bin/env perl

use File::Basename;
use File::Temp qw/ tempfile tempdir /;
use strict;
use Getopt::Long;

my $USAGE = "generateBareWigs.pl outputPrefix sample1REPLACECHROM.bam sample1REPLACECHROM.bam .. or ..\n" . 
"generateBareWigs.pl outputPrefix sample.bam .. or ..";
#GetOptions ('minct=i',\$minct, 'genome=s',\$genome, 'doCsv!' => \$doCsv, 'doOnlyRefCg!'=> \$doOnlyRefCg, 'doBare!' => \$doBare, 'intervalFile=s',\$intervalFile) || die "$USAGE ($doBare)\n";
#print STDERR "doBare=${doBare}\tintervalFile=${intervalFile}\n";

# Input params
die "$USAGE\n" unless (@ARGV==3);
my ($newname , $csvTemplate1, $csvTemplate2) = @ARGV;

# Constants
my $MATLAB = "/usr/usc/matlab/2011a/bin/matlab";
my $TAG = "GTF2MATLABMETHDIFF";
my $classpath = "/home/rcf-40/bberman/Java/uecgatk/bin:/home/rcf-40/bberman/Java/uecgatk/lib/GenomeAnalysisTK.jar:/home/rcf-40/bberman/Java/uecgatk/lib/genomeLibs.jar:/home/rcf-40/bberman/Java/uecgatk/lib/StingUtils.jar:/home/rcf-40/bberman/Java/uecgatk/lib/UscKeck.jar:/home/rcf-40/bberman/svn/genomeLibs/trunk/sam-1.26.jar:/home/rcf-40/bberman/svn/genomeLibs/trunk/charts4j-1.2.jar:/home/rcf-40/bberman/svn/genomeLibs/trunk/biojava-live_1.6/apps-live.jar:/home/rcf-40/bberman/svn/genomeLibs/trunk/biojava-live_1.6/biojava-live.jar:/home/rcf-40/bberman/svn/genomeLibs/trunk/biojava-live_1.6/bytecode.jar:/home/rcf-40/bberman/svn/genomeLibs/trunk/commons-math-1.1.jar:/home/rcf-40/bberman/svn/genomeLibs/trunk/biojava-live_1.6/commons-cli.jar:/home/rcf-40/bberman/svn/genomeLibs/trunk/biojava-live_1.6/commons-collections-2.1.jar:/home/rcf-40/bberman/svn/genomeLibs/trunk/biojava-live_1.6/commons-dbcp-1.1.jar:/home/rcf-40/bberman/svn/genomeLibs/trunk/biojava-live_1.6/commons-pool-1.1.jar:/home/rcf-40/bberman/svn/genomeLibs/trunk/biojava-live_1.6/demos-live.jar:/home/rcf-40/bberman/svn/genomeLibs/trunk/biojava-live_1.6/jgrapht-jdk1.5.jar:/home/rcf-40/bberman/svn/genomeLibs/trunk/biojava-live_1.6/junit-4.4.jar:/home/uec-00/shared/production/software/picard/default/picard-tools-1.40.zip";
my $matlabpath = "/home/uec-00/bberman/MatlabLibs/usc_matlab:/home/uec-00/bberman/MatlabLibs/insitu_matlab_g5:/home/uec-00/bberman/MatlabLibs/insitu_matlab";
# my $genomeToRef = {"hg18" => "hg18_unmasked", "hg19" => "hg19_rCRSchrm" };
# my $genomeToCpgsDir = {"hg18" => "Cpgls *.s", "hg19" => "CpgsHg19" };
# my $genomeRef = $genomeToRef->{$genome};
# die "Can't find genome ref for $genome\n" unless ($genomeRef);
# # Constants
# my $refDir = "/home/uec-00/shared/production/genomes/${genomeRef}";
# my $refFn = "${refDir}/${genomeRef}.fa";


my @regions = ("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", 
	       "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", 
	       "chr19", "chr20", "chr21", "chr22","chrX");
#@regions = ("chr2", "chr3", "chr5", "chr13");

@regions = ("chr11","chr4");
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
#my @wigPrefixes = ();
#my @bedgraphPrefixes = ();
foreach my $region (@regions)
{
#     # 0-(nMinCsVariants) variable wig, #nMinCsVariants=bare wig
#     for (my $round = 0; $round <= $endround; $round++)
#     {
	my $csvFn1 = ${csvTemplate1};
	$csvFn1 =~ s/REPLACECHROM/${region}/g if ($sepChrs);
	my $csvFn2 = ${csvTemplate2};
	$csvFn2 =~ s/REPLACECHROM/${region}/g if ($sepChrs);
	my $prefix = "${newname}-${region}";
	my $cmd = "export CLASSPATH=${classpath}; export MATLABPATH=${matlabpath} ; ";

	$cmd .= "echo \"a=\'${csvFn1}\'; b=\'${csvFn2}\'; processWGBSdiffMeth(a,b,\'normal\',\'tumor\',\'${prefix}\',5,35,[],0,50);\" ";
	$cmd .= " | ${MATLAB} -nodisplay -nosplash -logfile ${prefix}-logger.txt";

	my @dependJobs = ($holdJobId);
    
	my $curjobid = runCmd(0,$cmd, "${TAG}_matlab_r", \@dependJobs);
	push(@matlabJobs, $curjobid);
}


# # Now go through and concatenate
# my @catJobs = ();
# my @catOuts = ();
# # 0-(nMinCsVariants) variable wig, #nMinCsVariants=bare wig
# for (my $round = 0; $round <= $endround; $round++)
# {
#     print STDERR "\n\nOn Round $round\n";
#     my $mc = 0;
#     my $suffix;
#     my @extensions;
#     if ($round<$nMinCsVariants)
#     {
# 	$mc = @MIN_CS_VARIANTS[$round];
# 	$suffix = ($mc<100) ? "-CG-minc${mc}-maxw10000.variable" : "-CG-minwind${mc}.fixed";
# 	@extensions = ("wig", "bedGraph");
#     }
#     else
#     {
# 	$suffix = ".CG-minct${minct}-minconv1";
# 	@extensions = ($doCsv) ? ("csv") : ("wig");
#     }

#     foreach my $extension (@extensions)
#     {
# 	my ($cmd, $fnout) = concatCmd(\@regions, $newname, $suffix, $extension);
# 	my @dependJobs = @matlabJobs;
# 	my $curjobid = runCmd(0,$cmd, "${TAG}_concat_r${round}_e${extension}", \@dependJobs);
# 	push(@catJobs, $curjobid);
# 	push(@catOuts, $fnout);
#     }

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
my $nsecs = 1;
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
    if (0)
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

