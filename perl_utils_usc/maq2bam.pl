#!/usr/bin/env perl

use File::Basename;
use strict;


my $SAMDIR = "/home/uec-00/shared/production/software/samtools";

my $mapFn = $ARGV[0] || die "no input map";
my $refFa = $ARGV[1] || die "no input genome";

my $finalProcId = 0;
my $finalBamFn = 0;

die "input must be map file" if ($mapFn !~ /\.map$/);
die "not supported for non-human, will fix soon" if ($refFa !~ /hg/);
my $mapFnBase = basename($mapFn, qr/\.map/);

# Since we're going into a temp dir one level deeper, we have to adjust for relative FNs
my $curIn = $mapFn;
my $curOut;

# Maq 2 sam
$curOut = "${mapFnBase}.sam";
system "${SAMDIR}/maq2sam-long ${curIn} > ${curOut}";


# Sam to full bam
$curIn = $curOut;
$curOut = "${mapFnBase}.bam";
system "${SAMDIR}/samtools view -bt /home/uec-00/shared/production/genomes/sambam/hg18.fai -o ${curOut} ${curIn}";

# remove dups
$curIn = $curOut;
$curOut = "${mapFnBase}.NODUPS.bam";
system "${SAMDIR}/samtools rmdup -s ${curIn} ${curOut}";

# Sort bam
$curIn = $curOut;
$curOut = "${mapFnBase}.NODUPS.sorted.bam";
my $curOutPrefix = $curOut; $curOutPrefix =~ s/\.bam$//g; # You only specify the prefix
system "${SAMDIR}/samtools sort ${curIn} ${curOutPrefix}";

# calmd
$curIn = $curOut;
$curOut = "${mapFnBase}.NODUPS.sorted.calmd.bam";
system "${SAMDIR}/samtools calmd -b ${curIn} ${refFa} > ${curOut}";

# Set final BAM name
$finalBamFn = $curOut;

# Index bam
$curIn = $curOut;
$curOut = "${mapFnBase}.NODUPS.sorted.calmd.bam.bai";
system "${SAMDIR}/samtools index ${curIn} ${curOut}";
