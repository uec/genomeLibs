#!/usr/bin/env perl

use File::Basename;
use File::Temp qw/ tempfile tempdir /;
use strict;

my $USAGE = "generateMaqToBam.pl newname refGenome.fa f1.map f2.map f3.map ...";
my $region = "chr11";
my $SAMDIR = "/home/uec-00/bberman/bin";
my $RMTMPS = 1;
my $BATCHSIZE = 50;

# Input params
die "$USAGE\n" unless (@ARGV>=3);
my ($newname , $refFa, @mapFns) = @ARGV;

# Go through pipeline for each map.  We do it in batches because
# Map merge can only take a limited number.
my @topDependJobs = ();
my @topBamFns = ();
my @dependJobs = ();
my @bamFns = ();
my $onJob = 0;
my $tmpDir;
my $batchname = "";
foreach my $mapFn (@mapFns)
{
    if (($onJob % $BATCHSIZE) == 0)
    {
	# End old batch
	if (scalar(@dependJobs)>1)
	{
	    my ($mergeJobid, $mergeBamFn) = mergeBams($batchname, \@dependJobs, \@bamFns);
	    push(@topDependJobs, $mergeJobid);
	    push(@topBamFns, $mergeBamFn);
	}

	# Start a new batch
	@dependJobs = ();
	@bamFns = ();
	my $batchEnd = ($onJob+$BATCHSIZE-1);
	$batchname = "${newname}-BATCH${onJob}-${batchEnd}";
	$tmpDir = "./tmp_${batchname}";
	`mkdir $tmpDir`;
    }

    my ($jobid, $bamFn) = runMapPipeline($mapFn, $refFa, $tmpDir);
    push (@dependJobs, $jobid);
    push (@bamFns, "${tmpDir}/$bamFn");

    $onJob++;
}

# Finish final batch
if (scalar(@dependJobs)>1)
{
    my ($mergeJobid, $mergeBamFn) = mergeBams($batchname, \@dependJobs, \@bamFns);
    push(@topDependJobs, $mergeJobid);
    push(@topBamFns, $mergeBamFn);
}

# And merge top level
mergeBams($newname, \@topDependJobs, \@topBamFns);


# - - - - - - OUTPUTS


# Returns
# ($processId, $outBamFn)
sub mergeBams
{
    my ($prefix, $dependJobs, $individualBamFns) = @_;

    print STDERR "merging bams (".join(",",@$individualBamFns) .")\n\tFrom jobs(".join(",",@$dependJobs).")\n";

#    my $regionSec = ($region) ? ".${region}" : "";
    my $regionSec = "";
    my $curOut = "${prefix}.NODUPS.sorted.calmd${regionSec}.bam";
    my $curJobids;
    if (scalar(@$individualBamFns) == 1)
    {
	my $cmd = "mv " . @{$individualBamFns}[0] . " $curOut";
	$curJobids = [runCmd(0, $cmd, "M2B_mvBams", $dependJobs)];
    }
    else
    {
	my $cmd = "${SAMDIR}/samtools merge ${curOut} " . join(" ", @$individualBamFns);
	$cmd .= "; rm -f " . join(" ", @$individualBamFns) if ($RMTMPS);
	$curJobids = [runCmd(0, $cmd, "M2B_mergeBams", $dependJobs)];
    }

    # We have to remove dups again
    my $curIn = $curOut;
    $curOut = "${prefix}.NODUPS.sorted.calmd${regionSec}.NODUPS.bam";
    my $cmd = "${SAMDIR}/samtools rmdup -s ${curIn} ${curOut}";
    $cmd .= "; rm -f ${curIn}" if ($RMTMPS);
    $curJobids = [runCmd(0,$cmd, "M2B_rmdups", $curJobids)];

    my $finalBamOut = $curOut;

    # Index bam
    my $curIn = $curOut;
    $curOut = "${curIn}.bai";
    my $cmd = "${SAMDIR}/samtools index ${curIn} ${curOut}";
    $curJobids = [runCmd(0, $cmd, "M2B_index", $curJobids)];


    # Pull region
    if ($region)
    {
	$curIn = $finalBamOut;
	$curOut = "${prefix}.NODUPS.sorted.calmd.${region}.NODUPS.bam";
	my $cmd = "${SAMDIR}/samtools view -b -o ${curOut} ${curIn} ${region}";
	$curJobids = [runCmd(0,$cmd, "M2B_pullRegion", $curJobids)];

	# Don't do this. Make the region a side effect and actually merge
	# full bam going up
	# Set final BAM name
	#$finalBamFn = $curOut;

	# Index bam
	$curIn = $curOut;
	$curOut = "${curIn}.bai";
	my $cmd = "${SAMDIR}/samtools index ${curIn} ${curOut}";
	$curJobids = [runCmd(0,$cmd, "M2B_pullRegionIndex", $curJobids)];
    }

    return (@${curJobids}[0], $finalBamOut);
}


# Returns the final process ID and the BAM outfile
# ($processId, $outBamFn)
sub runMapPipeline
{
    my ($mapFn, $refFa, $tmpdir) = @_;

    my $finalProcId = 0;
    my $finalBamFn = 0;

    my $mapFnBase = basename($mapFn, qr/\.map/);

    # Since we're going into a temp dir one level deeper, we have to adjust for relative FNs
    my $curIn = $mapFn;
    if ($curIn =~ /^\.\//)
    {
	$curIn =~ s/^\.\//\.\.\//g;
    }
    elsif ($curIn =~ /^\//)
    {
	# Global, no problem
    }
    else
    {
	$curIn = "../" . $curIn;
    }


    my $curOut;
    my $curJobids = 0;

    # Maq 2 sam
    $curOut = "${mapFnBase}.sam";
    my $cmd = "${SAMDIR}/maq2sam-long ${curIn} > ${curOut}";
    $curJobids = [runCmd($tmpdir,$cmd, "M2B_maq2sam", $curJobids)];
    
    # Sam to full bam
    $curIn = $curOut;
    $curOut = "${mapFnBase}.bam";
    my $cmd = "${SAMDIR}/samtools view -bt /home/uec-00/shared/production/genomes/sambam/hg18.fai -o ${curOut} ${curIn}";
    $cmd .= "; rm -f ${curIn}" if ($RMTMPS);
    $curJobids = [runCmd($tmpdir,$cmd, "M2B_sam2fullbam", $curJobids)];

    # remove dups
    $curIn = $curOut;
    $curOut = "${mapFnBase}.NODUPS.bam";
    my $cmd = "${SAMDIR}/samtools rmdup -s ${curIn} ${curOut}";
    $cmd .= "; rm -f ${curIn}" if ($RMTMPS);
    $curJobids = [runCmd($tmpdir,$cmd, "M2B_rmdups", $curJobids)];
    
    # Sort bam
    $curIn = $curOut;
    $curOut = "${mapFnBase}.NODUPS.sorted.bam";
    my $curOutPrefix = $curOut; $curOutPrefix =~ s/\.bam$//g; # You only specify the prefix
    my $cmd = "${SAMDIR}/samtools sort ${curIn} ${curOutPrefix}";
    $cmd .= "; rm -f ${curIn}" if ($RMTMPS);
    $curJobids = [runCmd($tmpdir,$cmd, "M2B_sort", $curJobids)];

    # calmd
    $curIn = $curOut;
    $curOut = "${mapFnBase}.NODUPS.sorted.calmd.bam";
    my $cmd = "${SAMDIR}/samtools calmd -b ${curIn} ${refFa} > ${curOut}";
    $cmd .= "; rm -f ${curIn}" if ($RMTMPS);
    $curJobids = [runCmd($tmpdir,$cmd, "M2B_calmd", $curJobids)];

    # Set final BAM name
    $finalBamFn = $curOut;
    
    # Index bam
    $curIn = $curOut;
    $curOut = "${mapFnBase}.NODUPS.sorted.calmd.bam.bai";
    my $cmd = "${SAMDIR}/samtools index ${curIn} ${curOut}";
    $curJobids = [runCmd($tmpdir,$cmd, "M2B_index", $curJobids)];

    # Pull region
    if ($region)
    {
	$curIn = $finalBamFn;
	$curOut = "${mapFnBase}.NODUPS.sorted.calmd.${region}.bam";
	my $cmd = "${SAMDIR}/samtools view -b -o ${curOut} ${curIn} ${region}";
	$curJobids = [runCmd($tmpdir,$cmd, "M2B_pullRegion", $curJobids)];

	# Don't do this. Make the region a side effect and actually merge
	# full bam going up
	# # Set final BAM name
	# $finalBamFn = $curOut;

	# Index bam
	$curIn = $curOut;
	$curOut = "${mapFnBase}.NODUPS.sorted.calmd.${region}.bam.bai";
	my $cmd = "${SAMDIR}/samtools index ${curIn} ${curOut}";
	$curJobids = [runCmd($tmpdir,$cmd, "M2B_pullRegionIndex", $curJobids)];
    }

    $finalProcId = @${curJobids}[0];
    return ($finalProcId, $finalBamFn);
}


# Returns jobid
sub runCmd
{
    my ($tmpdir, $cmd, $prefix, $dependJobs) = @_;

    my ($fh, $file) = tempfile( "${prefix}XXXXXX" , DIR => "/tmp");
    print $fh "#Run on 1 processors on laird\n";
    print $fh "#PBS -l walltime=30:00:00\n#PBS -l mem=4000mb\n#PBS -l arch=x86_64\n#PBS -q laird\n";
    print $fh "#PBS -W depend=afterany:" . join(":",@$dependJobs) ."\n" if ($dependJobs && @$dependJobs);
    print $fh "cd \"\$PBS_O_WORKDIR\"\n";
    print $fh "${cmd}\n\n";
    close($fh);

    my $fullcmd = "qsub $file";
    $fullcmd = "cd ${tmpdir}; ${fullcmd}; cd .." if ($tmpdir);
    print STDERR "${cmd}\n${fullcmd}\n";

    my $lineout = `$fullcmd`;
    chomp $lineout;
    print STDERR "runCmd got lineout:\t\"$lineout\"\n";

    `echo "job ${lineout} =" >> M2B_pbs_log.txt`;
    `cat $file >> M2B_pbs_log.txt`;



    unlink($file);
    return $lineout;
}

