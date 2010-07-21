#!/usr/bin/env perl

use File::Basename;
use strict;


my $SAMDIR = "/home/uec-00/shared/production/software/samtools";

my $samFn = $ARGV[0] || die "no input sam file";
my $refFa = $ARGV[-1] || die "no input genome";
my $refFai = $refFa;
if( -e $refFa && $refFa =~ /\.fa$/)
{
                $refFai = $refFai . ".fai";
                die "no fai found" unless -e $refFai;
}
else
{
        die "no ref or ref index found";
}

#handle multiple sams for a merge
if(scalar(@ARGV) > 2)
{
        my @checkedSamFiles;
        for my $samfile (@ARGV[0..($#ARGV - 1)])
        {
                if(-e $samfile)
                {
                        push @checkedSamFiles, $samfile;
                }
                else
                {
                        print STDERR "$samfile NOT FOUND!";
                } 
        }
        my $mergeFileString = join(" ", @checkedSamFiles);
        my $mergedFileName = $samFn;
        $mergedFileName = basename($mergedFileName);
        $mergedFileName =~ s/\.sam$/\.merged\.sam/;
        print STDERR "cat $mergeFileString > $mergedFileName\n";
        system("cat $mergeFileString > $mergedFileName");
        $samFn = $mergedFileName;
}


my $finalProcId = 0;
my $finalBamFn = 0;

die "input must be sam file" if ($samFn !~ /\.sam$/);
my $samFnBase = basename($samFn);
$samFnBase =~ s/\.sam$//;

# Since we're going into a temp dir one level deeper, we have to adjust for relative FNs
my $curIn = $samFn;



# Sam to full bam
my $curOut = "${samFnBase}.bam";
print STDERR "${SAMDIR}/samtools view -bt $refFai -o ${curOut} ${curIn}\n";
system "${SAMDIR}/samtools view -bt $refFai -o ${curOut} ${curIn}";
exit;

# remove dups
$curIn = $curOut;
$curOut = "${samFnBase}.NODUPS.bam";
print STDERR "${SAMDIR}/samtools rmdup -s ${curIn} ${curOut}\n";
system "${SAMDIR}/samtools rmdup -s ${curIn} ${curOut}";

# Sort bam
$curIn = $curOut;
$curOut = "${samFnBase}.NODUPS.sorted.bam";
my $curOutPrefix = $curOut; $curOutPrefix =~ s/\.bam$//g; # You only specify the prefix
print STDERR "${SAMDIR}/samtools sort ${curIn} ${curOutPrefix}\n";
system "${SAMDIR}/samtools sort ${curIn} ${curOutPrefix}";

# calmd
$curIn = $curOut;
$curOut = "${samFnBase}.NODUPS.sorted.calmd.bam";
print STDERR "${SAMDIR}/samtools calmd -b ${curIn} ${refFa} > ${curOut}\n";
system "${SAMDIR}/samtools calmd -b ${curIn} ${refFa} > ${curOut}";

# Set final BAM name
$finalBamFn = $curOut;

# Index bam
$curIn = $curOut;
$curOut = "${samFnBase}.NODUPS.sorted.calmd.bam.bai";
print STDERR "${SAMDIR}/samtools index ${curIn} ${curOut}\n";
system "${SAMDIR}/samtools index ${curIn} ${curOut}";