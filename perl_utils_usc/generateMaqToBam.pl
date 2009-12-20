#!/usr/bin/env perl

use File::Basename;
use File::Temp qw/ tempfile tempdir /;
use strict;

my $USAGE = "generateMaqToBam.pl newname refGenome.fa f1.map f2.map f3.map ...";
my $region = "chr11";
my $SAMDIR = "/home/uec-00/bberman/bin";
my $RMTMPS = 1;

# Input params
die "$USAGE\n" unless (@ARGV>=3);
my ($newname , $refFa, @mapFns) = @ARGV;


# Go through pipeline for each map
my @dependJobs = ();
my @individualBamFns = ();
foreach my $mapFn (@mapFns)
{
    my ($jobid, $bamFn) = runMapPipeline($mapFn, $refFa);
    push (@dependJobs, $jobid);
    push (@individualBamFns, $bamFn);
}

mergeBams($newname, \@dependJobs, \@individualBamFns);


# - - - - - - OUTPUTS

sub mergeBams
{
    my ($prefix, $dependJobs, $individualBamFns) = @_;

    print STDERR "merging bams (".join(",",@$individualBamFns) .")\n\tFrom jobs(".join(",",@$dependJobs).")\n";

    my $regionSec = ($region) ? ".${region}" : "";
    my $curOut = "${prefix}.NODUPS.sorted.calmd${regionSec}.bam";
    my $cmd = "${SAMDIR}/samtools merge ${curOut} " . join(" ", @$individualBamFns);
    my $curJobids = [runCmd($cmd, "M2B_mergeBams", $dependJobs)];

    # Index bam
    my $curIn = $curOut;
    $curOut = "${curIn}.bai";
    my $cmd = "${SAMDIR}/samtools index ${curIn} ${curOut}";
    $curJobids = [runCmd($cmd, "M2B_index", $curJobids)];
}


# Returns the final process ID and the BAM outfile
# ($processId, $outBamFn)
sub runMapPipeline
{
    my ($mapFn, $refFa) = @_;

    my $finalProcId = 0;
    my $finalBamFn = 0;

    my $mapFnBase = basename($mapFn, qr/\.map/);

    my $curIn = $mapFn;
    my $curOut;
    my $curJobids = 0;

    `mkdir ./tmp`;

    # Maq 2 sam
    $curOut = "./tmp/${mapFnBase}.sam";
    my $cmd = "${SAMDIR}/maq2sam-long ${curIn} > ${curOut}";
    $curJobids = [runCmd($cmd, "M2B_maq2sam", $curJobids)];
    
    # Sam to full bam
    $curIn = $curOut;
    $curOut = "./tmp/${mapFnBase}.bam";
    my $cmd = "${SAMDIR}/samtools view -bt /home/uec-00/shared/production/genomes/sambam/hg18.fai -o ${curOut} ${curIn}";
    $cmd .= "; rm -f ${curIn}" if ($RMTMPS);
    $curJobids = [runCmd($cmd, "M2B_sam2fullbam", $curJobids)];

    # remove dups
    $curIn = $curOut;
    $curOut = "./tmp/${mapFnBase}.NODUPS.bam";
    my $cmd = "${SAMDIR}/samtools rmdup -s ${curIn} ${curOut}";
    $cmd .= "; rm -f ${curIn}" if ($RMTMPS);
    $curJobids = [runCmd($cmd, "M2B_rmdups", $curJobids)];
    
    # Sort bam
    $curIn = $curOut;
    $curOut = "./tmp/${mapFnBase}.NODUPS.sorted.bam";
    my $curOutPrefix = $curOut; $curOutPrefix =~ s/\.bam$//g; # You only specify the prefix
    my $cmd = "${SAMDIR}/samtools sort ${curIn} ${curOutPrefix}";
    $cmd .= "; rm -f ${curIn}" if ($RMTMPS);
    $curJobids = [runCmd($cmd, "M2B_sort", $curJobids)];

    # calmd
    $curIn = $curOut;
    $curOut = "./tmp/${mapFnBase}.NODUPS.sorted.calmd.bam";
    my $cmd = "${SAMDIR}/samtools calmd -b ${curIn} ${refFa} > ${curOut}";
    $cmd .= "; rm -f ${curIn}" if ($RMTMPS);
    $curJobids = [runCmd($cmd, "M2B_calmd", $curJobids)];

    # Set final BAM name
    $finalBamFn = $curOut;
    
    # Index bam
    $curIn = $curOut;
    $curOut = "./tmp/${mapFnBase}.NODUPS.sorted.calmd.bam.bai";
    my $cmd = "${SAMDIR}/samtools index ${curIn} ${curOut}";
    $curJobids = [runCmd($cmd, "M2B_index", $curJobids)];

    # Pull region
    if ($region)
    {
	$curIn = "./tmp/${mapFnBase}.NODUPS.sorted.calmd.bam";
	$curOut = "./tmp/${mapFnBase}.NODUPS.sorted.calmd.${region}.bam";
	my $cmd = "${SAMDIR}/samtools view -b -o ${curOut} ${curIn} ${region}";
	$curJobids = [runCmd($cmd, "M2B_pullRegion", $curJobids)];

	# Set final BAM name
	$finalBamFn = $curOut;

	# Index bam
	$curIn = $curOut;
	$curOut = "./tmp/${mapFnBase}.NODUPS.sorted.calmd.${region}.bam.bai";
	my $cmd = "${SAMDIR}/samtools index ${curIn} ${curOut}";
	$curJobids = [runCmd($cmd, "M2B_pullRegionIndex", $curJobids)];
    }

    $finalProcId = @${curJobids}[0];
    return ($finalProcId, $finalBamFn);
}


# Returns jobid
sub runCmd
{
    my ($cmd, $prefix, $dependJobs) = @_;

    my ($fh, $file) = tempfile( "${prefix}XXXXXX" , DIR => "/tmp");
    print $fh "#Run on 1 processors on laird\n";
    print $fh "#PBS -l walltime=30:00:00\n#PBS -l mem=2000mb\n#PBS -l arch=x86_64\n#PBS -q laird\n";
    print $fh "#PBS -W depend=afterany:" . join(":",@$dependJobs) ."\n" if ($dependJobs && @$dependJobs);
    print $fh "cd \"\$PBS_O_WORKDIR\"\n";
    print $fh "${cmd}\n\n";
    close($fh);


    my $fullcmd = "qsub $file";
    print STDERR "${cmd}\n${fullcmd}\n";
    my $lineout = `$fullcmd`;
    chomp $lineout;
    print STDERR "runCmd got lineout:\t\"$lineout\"\n";

    `echo "job ${lineout} =" >> M2B_pbs_log.txt`;
    `cat $file >> M2B_pbs_log.txt`;



    unlink($file);
    return $lineout;
}

