#!/usr/bin/env perl

use File::Basename;
use strict;


my $SAMDIR = "/home/uec-00/shared/production/software/samtools";
my $PICARDDIR = "/home/uec-00/shared/production/software/picard/default/";
my $JAVA = "/home/uec-00/shared/production/software/java/default/bin/java";

my $mapFn = $ARGV[0] || die "no input map";
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

die "input must be map file" if ($mapFn !~ /\.map$/);
#die "not supported for non-human, will fix soon" if ($refFa !~ /hg/);
my $mapFnBase = basename($mapFn);
$mapFnBase =~ s/\.map$//;

# Since we're going into a temp dir one level deeper, we have to adjust for relative FNs
my $curIn = $mapFn;
my $curOut;

# Maq 2 sam
$curOut = "${mapFnBase}.sam";
system "${SAMDIR}/maq2sam-long ${curIn} > ${curOut}";


# Sam to full bam
$curIn = $curOut;
$curOut = "${mapFnBase}.bam";
system "${SAMDIR}/samtools view -bt $refFai -o ${curOut} ${curIn}";


# remove dups (legacy mode for bw compat)
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
system "${SAMDIR}/samtools calmd -b ${curIn} ${refFa} > ${curOut} 2> /dev/null";

# Set final BAM name
$finalBamFn = $curOut;

# Index bam
$curIn = $curOut;
$curOut = "${mapFnBase}.NODUPS.sorted.calmd.bam.bai";
system "${SAMDIR}/samtools index ${curIn} ${curOut}";

# sort with picard (no point on making a nodups file seperately, but do it for now)
$curIn = "${mapFnBase}.bam";
$curOut = "${mapFnBase}.picard.sorted.bam";
system "$JAVA -Xmx14g -jar $PICARDDIR/SortSam.jar INPUT=$curIn OUTPUT=$curOut SORT_ORDER=coordinate MAX_RECORDS_IN_RAM=3000000";

# mark dups with picard (no point on making a nodups file seperately, but do it for now)
$curIn = $curOut;
$curOut = "${mapFnBase}.picard.nodups.bam";
system "$JAVA -Xmx14g -jar $PICARDDIR/MarkDuplicates.jar INPUT=$curIn OUTPUT=$curOut METRICS_FILE=met.txt MAX_RECORDS_IN_RAM=3000000";

#replace original bam with dup-marked, sorted, bam
$curIn = $curOut;
$curOut = "${mapFnBase}.bam";
system "mv $curIn $curOut";
