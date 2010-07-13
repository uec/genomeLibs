#!/usr/bin/env perl

use File::Basename;
use strict;


my $SAMDIR = "/home/uec-00/shared/production/software/samtools";

my $samFn = $ARGV[0] || die "no input sam file";
my $refFa = $ARGV[1] || die "no input genome";
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

my $finalProcId = 0;
my $finalBamFn = 0;

die "input must be map file" if ($samFn !~ /\.sam$/);
my $samFnBase = basename($samFn);
$samFnBase =~ s/\.sam$//;

# Since we're going into a temp dir one level deeper, we have to adjust for relative FNs
my $curIn = $samFn;
my $curOut;



# Sam to full bam
$curIn = $curOut;
$curOut = "${samFnBase}.bam";
system "${SAMDIR}/samtools view -bt $refFai -o ${curOut} ${curIn}";

# remove dups
$curIn = $curOut;
$curOut = "${samFnBase}.NODUPS.bam";
system "${SAMDIR}/samtools rmdup -s ${curIn} ${curOut}";

# Sort bam
$curIn = $curOut;
$curOut = "${samFnBase}.NODUPS.sorted.bam";
my $curOutPrefix = $curOut; $curOutPrefix =~ s/\.bam$//g; # You only specify the prefix
system "${SAMDIR}/samtools sort ${curIn} ${curOutPrefix}";

# calmd
$curIn = $curOut;
$curOut = "${samFnBase}.NODUPS.sorted.calmd.bam";
system "${SAMDIR}/samtools calmd -b ${curIn} ${refFa} > ${curOut}";

# Set final BAM name
$finalBamFn = $curOut;

# Index bam
$curIn = $curOut;
$curOut = "${samFnBase}.NODUPS.sorted.calmd.bam.bai";
system "${SAMDIR}/samtools index ${curIn} ${curOut}";
