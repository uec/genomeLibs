#!/usr/bin/env perl

use File::Basename;
use File::Temp qw/ tempfile tempdir /;
use strict;
use Getopt::Long;
use File::Basename;
use lib dirname (__FILE__);
use EpigenomeUtils;

my $doBare = 0;
my $intervalFile = 0; # If we pass these in , they're done in parallel
my $doCsv = 0; # For
my $doOnlyRefCg = 0; # Speeds it up considerably
my $genome = "hg18";
my $minct = 1;
my $USAGE = "generateBareWigs.pl minCpgs|0 outputPrefix sampleREPLACECHROM.bam .. or ..\n" . 
"generateBareWigs.pl minCpgs outputPrefix sample.bam .. or ..";
GetOptions ('minct=i',\$minct, 'genome=s',\$genome, 'doCsv!' => \$doCsv, 'doOnlyRefCg!'=> \$doOnlyRefCg, 'doBare!' => \$doBare, 'intervalFile=s',\$intervalFile) || die "$USAGE ($doBare)\n";
print STDERR "doBare=${doBare}\tintervalFile=${intervalFile}\n";

# Input params
die "$USAGE\n" unless (@ARGV==3);
my ($mincs, $newname , $bamTemplate) = @ARGV;
my @MIN_CS_VARIANTS = ($mincs>0) ? ($mincs) : ();

# Constants
my $genomeToRef = {"hg18" => "hg18_unmasked", "hg19" => "hg19_rCRSchrm" };
my $genomeToCpgsDir = {"hg18" => "Cpgs", "hg19" => "CpgsHg19" };
my $genomeRef = $genomeToRef->{$genome};
die "Can't find genome ref for $genome\n" unless ($genomeRef);

# Constants
my $TAG = "BAM2WIG";
my $IGVTOOLS = "$SOFTWAREROOT/software/igvtools/igvtools";
my $CGINTERVALS = "$GENOMEROOT/genomic-data-misc/CpgsAll.${genome}.bed";
my $CGINTERVALS_SEP = "$GENOMEROOT/genomic-data-misc/Cpgs/${genomeRef}/REPLACECHROM.fa.cpgCoords.bed";
my $refDir = "$GENOMEROOT/genomes/${genomeRef}";
my $refFn = "${refDir}/${genomeRef}.fa";


#my @SUFFIXES = qw/GCH HCG readcvg/;


my @regions = ("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", 
	       "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", 
	       "chr19", "chr20", "chr21", "chr22"); #, "chrX");
#@regions = ("chr2", "chr3", "chr5", "chr13");

@regions = ("chr11","chr4");
#my $CGINTERVALS_SEP = "/home/uec-02/bberman/BSseq/tumor/genomic-data-misc/Cpgs/${genomeRef}/CpgsRef.REPLACECHROM.7M-9M.bed";

my $classpath = "/home/rcf-40/bberman/Java/uecgatk/bin:/home/rcf-40/bberman/Java/uecgatk/lib/GenomeAnalysisTK.jar:/home/rcf-40/bberman/Java/uecgatk/lib/genomeLibs.jar:/home/rcf-40/bberman/Java/uecgatk/lib/StingUtils.jar:/home/rcf-40/bberman/Java/uecgatk/lib/UscKeck.jar:/home/rcf-40/bberman/svn/genomeLibs/trunk/sam-1.26.jar:/home/rcf-40/bberman/svn/genomeLibs/trunk/charts4j-1.2.jar:/home/rcf-40/bberman/svn/genomeLibs/trunk/biojava-live_1.6/apps-live.jar:/home/rcf-40/bberman/svn/genomeLibs/trunk/biojava-live_1.6/biojava-live.jar:/home/rcf-40/bberman/svn/genomeLibs/trunk/biojava-live_1.6/bytecode.jar:/home/rcf-40/bberman/svn/genomeLibs/trunk/commons-math-1.1.jar:/home/rcf-40/bberman/svn/genomeLibs/trunk/biojava-live_1.6/commons-cli.jar:/home/rcf-40/bberman/svn/genomeLibs/trunk/biojava-live_1.6/commons-collections-2.1.jar:/home/rcf-40/bberman/svn/genomeLibs/trunk/biojava-live_1.6/commons-dbcp-1.1.jar:/home/rcf-40/bberman/svn/genomeLibs/trunk/biojava-live_1.6/commons-pool-1.1.jar:/home/rcf-40/bberman/svn/genomeLibs/trunk/biojava-live_1.6/demos-live.jar:/home/rcf-40/bberman/svn/genomeLibs/trunk/biojava-live_1.6/jgrapht-jdk1.5.jar:/home/rcf-40/bberman/svn/genomeLibs/trunk/biojava-live_1.6/junit-4.4.jar:/home/uec-00/shared/production/software/picard/default/picard-tools-1.40.zip";

#my $minCphFrac = 1.01;

my $lastChr = "";
my $lastEnd = 0;
if ($intervalFile)
{
    die "Couldn't read $intervalFile\n" unless open(INTF, $intervalFile);
    @regions = ();
    while (my $line = <INTF>)
    {
	chomp $line;
	next if ($line =~ /^\s*\#/);
	next if ($line =~ /^\s*$/);
	my ($chr, $s, $e, $strand, $name) = split(/\t/,$line);
	$s = $s-200000;
	$s = 1 if ($s<1);
	$e = ($e+200000);
	die "Overlapping interval $name\n" if ( ($chr eq $lastChr) && ($s <= $lastEnd));
	my $reg = sprintf("%s:%d-%d",$chr,$s,$e);
	push(@regions, $reg);
	$lastChr = $chr;
	$lastEnd = $e;
	print STDERR "$reg\n";
    }
    close(INTF);
}


#my @MIN_CS_VARIANTS = (10, 50, 100);
my $nMinCsVariants = scalar(@MIN_CS_VARIANTS);

# special case.  If newname contains a full path filename, remove it
$newname = basename($newname);

# Go through and check files
my $sepChrBams = ($bamTemplate =~ /REPLACECHROM/);
foreach my $region (@regions)
{
    my $fn = $bamTemplate;
    $fn =~ s/REPLACECHROM/${region}/g if ($sepChrBams);
    die "Following file does not exist. Quitting\n${fn}\n" unless (-f $fn);

}

# Start a fake command just to hold all jobs until the end
my $holdJobId = runCmd(0, "ls", "${TAG}_HOLDER", [], 1);
print STDERR "HOlding all jobs on job $holdJobId\n";


my $endround = ($doBare) ? ($nMinCsVariants) : ($nMinCsVariants-1);
my @wigJobs = ();
#my @wigPrefixes = ();
#my @bedgraphPrefixes = ();
foreach my $region (@regions)
{
    # 0-(nMinCsVariants) variable wig, #nMinCsVariants=bare wig
    for (my $round = 0; $round <= $endround; $round++)
    {
	print STDERR "\n\nOn Round $round/$nMinCsVariants\n";
	my $mc = 0;
	$mc = @MIN_CS_VARIANTS[$round] if ($round<$nMinCsVariants);

	my $bamFn = ${bamTemplate};
	$bamFn =~ s/REPLACECHROM/${region}/g if ($sepChrBams);
	my $prefix = "${newname}-${region}";
	my $cmd = "export CLASSPATH=${classpath}";

	my $intervalsSec = "--intervals $region"; # -L region because we may use whole-genome BAM
	if ($doOnlyRefCg)
	{
	    my $sepRegions = ($region =~ /^chr/);
	    my $cgfile = ($sepRegions) ? $CGINTERVALS_SEP : $CGINTERVALS;
	    $cgfile =~ s/REPLACECHROM/${region}/g if ($sepRegions);
	    $intervalsSec = "--intervals $cgfile";
	}
	my $csvSec = ($doCsv) ? "--csvMode" : "";
	

	if ($round<$nMinCsVariants)
	{
	    # If it's over 100 , treat it as a fixed wind
	    if ($mc > 100)
	    {
		$cmd .= " ; java -Xmx7500m org.broadinstitute.sting.gatk.CommandLineGATK -T BisulfiteSeqToCytosineVariableWig -R ${refFn} -I ${bamFn} -nt 1 --minCT ${minct} --minConv 1 -maxa 0.101 --minCpgs 10 --useFixedWind --fixedWindMinLength ${mc} --outPrefix ${prefix} --min_mapping_quality_score 30 -et NO_ET ${intervalsSec} ${csvSec}";
	    }
	    else
	    {
		$cmd .= " ; java -Xmx7500m org.broadinstitute.sting.gatk.CommandLineGATK -T BisulfiteSeqToCytosineVariableWig -R ${refFn} -I ${bamFn} -nt 1 --minCT ${minct} --minConv 1 -maxa 0.101 --minCpgs ${mc} --maxWindStretch 10000 --outPrefix ${prefix} --min_mapping_quality_score 30 -et NO_ET ${intervalsSec} ${csvSec}";
	    }
	}
	else
	{
	    # Do one wig and one csv
	    $cmd .= " ; java -Xmx7500m org.broadinstitute.sting.gatk.CommandLineGATK -T BisulfiteSeqToBareWig -R ${refFn} -I ${bamFn} -nt 1 --minCT ${minct} --minConv 1 -maxa 0.101 --outPrefix ${prefix} --min_mapping_quality_score 30 -et NO_ET ${intervalsSec}";
	    $cmd .= " ; java -Xmx7500m org.broadinstitute.sting.gatk.CommandLineGATK -T BisulfiteSeqToBareWig -R ${refFn} -I ${bamFn} -nt 1 --minCT ${minct} --minConv 1 -maxa 0.101 --outPrefix ${prefix} --min_mapping_quality_score 30 -et NO_ET ${intervalsSec} -doCsv";
	}
	
	my @dependJobs = ($holdJobId);
    
	my $curjobid = runCmd(0,$cmd, "${TAG}_wiggatk_r${round}", \@dependJobs);
	push(@wigJobs, $curjobid);
#	push(@wigPrefixes, $prefix);
#	push(@bedgraphPrefixes, $prefix) if ($round<$nMinCsVariants);
    }
}

if (0)
{
# Now go through and concatenate
my @catJobs = ();
my @catOuts = ();
# 0-(nMinCsVariants) variable wig, #nMinCsVariants=bare wig
for (my $round = 0; $round <= $endround; $round++)
{
    print STDERR "\n\nOn Round $round\n";
    my $mc = 0;
    my $suffix;
    my @extensions;
    if ($round<$nMinCsVariants)
    {
	$mc = @MIN_CS_VARIANTS[$round];
	$suffix = ($mc<100) ? "-CG-minc${mc}-maxw10000.variable" : "-CG-minwind${mc}.fixed";
	@extensions = ("wig", "bedGraph");
    }
    else
    {
	$suffix = ".CG-minct${minct}-minconv1";
	@extensions = ($doCsv) ? ("csv") : ("wig");
    }

    foreach my $extension (@extensions)
    {
	my ($cmd, $fnout) = concatCmd(\@regions, $newname, $suffix, $extension);
	my @dependJobs = @wigJobs;
	my $curjobid = runCmd(0,$cmd, "${TAG}_concat_r${round}_e${extension}", \@dependJobs);
	push(@catJobs, $curjobid);
	push(@catOuts, $fnout);
    }
}

# Now igvtools and gzip
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

