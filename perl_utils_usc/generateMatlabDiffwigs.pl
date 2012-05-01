#!/usr/bin/env perl

use File::Basename;
use File::Temp qw/ tempfile tempdir /;
use strict;
use Getopt::Long;

my $USAGE = "generateMatlabDiffwigs.pl outputPrefix numCompareSamples compareSample1REPLACECHROM.csv compareSample1REPLACECHROM.csv testSample1REPLACECHROM.csv testSample2REPLACECHROM.csv";
my $genome = "hg19";
my $compareDirection = 1;  # -1 for max, +1 for min
my $compareCutoff = 50;
my $testCutoff = -1;
my $testDirection = 1;  # -1 for max, +1 for min
my $compMinCount = 2;
my $testMinCount = 5;
my $resolution = 0;
GetOptions ('compareCutoff=f',\$compareCutoff,'testCutoff=f',\$testCutoff, 
	    'compareDirection=i',\$compareDirection, 'testDirection=i',\$testDirection, 
	    'genome=s',\$genome, 'resolution=i',\$resolution,
	    'compMinCount=i',\$compMinCount, 'testMinCount=i',\$testMinCount) || die "$USAGE\n";
#print STDERR "doBare=${doBare}\tintervalFile=${intervalFile}\n";

# Input params
die "$USAGE\n" unless (@ARGV>=3);
my ($newname , $numCompSamples, @csvTemplates) = @ARGV;

my @compCsvs = @csvTemplates[0..($numCompSamples-1)];
my @testCsvs = @csvTemplates[$numCompSamples..(scalar(@csvTemplates)-1)];

# Constants
my $MATLAB = "/usr/usc/matlab/2011a/bin/matlab";
my $IGVTOOLS = "/home/uec-00/shared/production/software/igvtools/default/igvtools";
my $TAG = "GTF2MATLABDIFFWIGS";
my $GL = "/home/rcf-40/bberman/svn/genomeLibs/trunk";
my $classpath = "/home/rcf-40/bberman/Java/uecgatk/bin:/home/rcf-40/bberman/Java/uecgatk/lib/GenomeAnalysisTK.jar:/home/rcf-40/bberman/Java/uecgatk/lib/genomeLibs.jar:/home/rcf-40/bberman/Java/uecgatk/lib/StingUtils.jar:/home/rcf-40/bberman/Java/uecgatk/lib/UscKeck.jar:${GL}/sam-1.26.jar:${GL}/charts4j-1.2.jar:${GL}/biojava-live_1.6/apps-live.jar:${GL}/biojava-live_1.6/biojava-live.jar:${GL}/biojava-live_1.6/bytecode.jar:${GL}/commons-math-1.1.jar:${GL}/biojava-live_1.6/commons-cli.jar:${GL}/biojava-live_1.6/commons-collections-2.1.jar:${GL}/biojava-live_1.6/commons-dbcp-1.1.jar:${GL}/biojava-live_1.6/commons-pool-1.1.jar:${GL}/biojava-live_1.6/demos-live.jar:${GL}/biojava-live_1.6/jgrapht-jdk1.5.jar:${GL}/biojava-live_1.6/junit-4.4.jar:/home/uec-00/shared/production/software/picard/default/picard-tools-1.40.zip:${GL}/lib/tuple.jar";
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
	       "chr19", "chr20", "chr21", "chr22", "chrX");
#@regions = ("chr2", "chr3", "chr5", "chr13");
@regions = ("chr8","chr21","chr22");


my $lastChr = "";
my $lastEnd = 0;

# special case.  If newname contains a full path filename, remove it
$newname = basename($newname);

# Go through and check files


my $sepChrs = ($csvTemplates[0] =~ /REPLACECHROM/);
foreach my $region (@regions)
{
    for my $csvTemplate (@csvTemplates)
    {
	my $fn = $csvTemplate;
	$fn =~ s/REPLACECHROM/${region}/g if ($sepChrs);
	die "Following file does not exist. Quitting\n${fn}\n" unless (-f $fn);
    }
}

# Start a fake command just to hold all jobs until the end
my $holdJobId = runCmd(0, "ls", "${TAG}_HOLDER", [], 1);
print STDERR "HOlding all jobs on job $holdJobId\n";

my @dependJobs;

# combine CpGs
@dependJobs = ($holdJobId);
my @combineJobs = ();
foreach my $region (@regions)
{
    for my $csvTemplate (@csvTemplates)
    {
	my $fn = $csvTemplate;
	$fn =~ s/REPLACECHROM/${region}/g if ($sepChrs);
	die "Following file does not exist. Quitting\n${fn}\n" unless (-f $fn);

	my $cmd = "matlabCsvCombineCpgs.pl ${fn}";

	my $curjobid = runCmd(0,$cmd, "${TAG}_combineCpGs", \@dependJobs);
	push(@combineJobs, $curjobid);
    }
}
@csvTemplates = map {s/\.csv$/\.combinedCG\.csv/g; $_} @csvTemplates;


# Make a fake one
my $fakeJobId = runCmd(0, "sleep 2", "${TAG}_FAKE", \@combineJobs);


@dependJobs = ($fakeJobId);
my @matlabJobs = ();
#foreach my $testCsvInd (${numCompSamples}..(scalar(@csvTemplates)-1))
# Do all including the comp samples
foreach my $testCsvInd (0..(scalar(@csvTemplates)-1))
{
    print STDERR "On testCsvInd=${testCsvInd}\n";

    my $suffix = "";
    my $prefix = "";

foreach my $region (@regions)
{
    my @bases = ();
    my @ids = ();
    my @csvFns = ();
    for (my $i = 0; $i < scalar(@csvTemplates); $i++)
    {
	$csvFns[$i] = $csvTemplates[$i];
	$csvFns[$i] =~ s/REPLACECHROM/${region}/g if ($sepChrs);
	$bases[$i] = basename($csvFns[$i]);
	$ids[$i] = idFromTCGA($bases[$i]);
    }

# FILE NAMING
    $prefix = sprintf("%s-%s",$newname,$ids[$testCsvInd]);
    my $cmd = "";
    $cmd .= "export CLASSPATH=${classpath}; export MATLABPATH=${matlabpath} ; ";

    my $filterDesc = sprintf("%s%.2f-%s%.2f-compmin%d-testmin%d", ($compareDirection>0)?"min":"max",$compareCutoff,
			     ($testDirection>0)?"min":"max",$testCutoff,$compMinCount,$testMinCount);
    $suffix = "-${filterDesc}";
    my $mainti = sprintf("%s-%s%s",$prefix,$region,$suffix);

    # MATLAB PARAMS
    my $aSec = "{" . join(", ", map {"\'$_\'"} @csvFns[0..($numCompSamples-1)]) . "}";
    my $bSec = "\'".$csvFns[$testCsvInd]."\'";
    $cmd .= "echo \"a=${aSec}; b=${bSec}; ";
    my $anSec = "{" . join(", ", map {"\'$_\'"} @bases[0..($numCompSamples-1)]) . "}";
    my $bnSec = "\'".$bases[$testCsvInd]."\'";
    $cmd .= " an=${anSec}; bn=${bnSec}; ";
    $cmd .= " processWGBSmakeFilteredWigs(a,b,an,bn,\'${mainti}\',$compareCutoff,$testCutoff,[],$compMinCount,$testMinCount,$resolution,$compareDirection,$testDirection,0);\" ";
    $cmd .= " | ${MATLAB} -nodisplay -nosplash -logfile ${prefix}-logger.txt";
    
    print STDERR "$cmd\n";

    # It adds stuff on
    $suffix .= ".${resolution}bpResolution";
#    my $gtfBase = "${newname}-${region}${suffix}";
    
    my $curjobid = runCmd(0,$cmd, "${TAG}_matlab_r", \@dependJobs);
    push(@matlabJobs, $curjobid);
}



# Concatenate
 my @catJobs = ();
 my @catOuts = ();
if (scalar(@regions)>0)
{
 	my ($cmd, $fnout) = concatCmd(\@regions, $prefix, $suffix, "wig");
 	my @dependJobs = @matlabJobs;
 	my $curjobid = runCmd(0,$cmd, "${TAG}_concat", \@dependJobs);
 	push(@catJobs, $curjobid);
 	push(@catOuts, $fnout);
}

#Now igvtools and gzip
foreach my $wigfn (@catOuts)
{
    my @dependJobs = @catJobs;
    my $tdffn = $wigfn;
    my @gzipDepends = @dependJobs;
    if ($tdffn =~ /\.wig/)
    {
	$tdffn =~ s/\.wig/\.tdf/g;
	my $cmd = "${IGVTOOLS} tile -c -f mean,min,max,median ${wigfn} ${tdffn} ${genome}";
	my $curjobid = runCmd(0,$cmd, "${TAG}_igvtools", \@dependJobs);
	push(@gzipDepends, $curjobid);
    }

    # Don't start the gzip until the tdf is finished, otherwise uncompressed file will disappear
    my $cmd = "gzip ${wigfn}";
    my $curjobid = runCmd(0,$cmd, "${TAG}_gzip", \@gzipDepends);
}


}


# Now we're ready to release the hounds
my $nsecs = 5;
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
	my $grepsec = ($count == 0) ? "" : " | grep -v track | grep -v browser ";
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
    print $fh "#PBS -l walltime=30:00:00\n#PBS -l mem=15995mb\n#PBS -l arch=x86_64\n#PBS -q lairdprio\n";
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

sub idFromTCGA
{
    my ($fullName) = @_;

    my $out = $fullName;
    if ($out =~ /(TCGA.*)\-05_/)
    {
	$out = $1;
    }
    return $out;
}
