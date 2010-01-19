#!/usr/bin/env perl

use strict;
use File::Basename;
use Getopt::Long;
use File::Temp qw/ tempfile /;

my $USAGE = "encodePeakFilesToMethylStats.pl -centeredSize 500 peaks1.peaks peaks2.peaks ...";


my $centeredSize = 0;
GetOptions ('centeredSize=i'=>\$centeredSize)|| die "$USAGE\n";

die "$USAGE\n" unless (@ARGV>0);

foreach my $peakFn (@ARGV)
{
    my ($name, $path, $suf) = fileparse($peakFn, qr/\.txt/); #qr/\.[^.]*/);

	#broadPeakFilesToGtf.pl ./wgEncodeYaleChIPseqPeaksHct116Tcf4V2.narrowPeak
	#java GtfFilterNonoverlapping ./wgEncodeYaleChIPseqPeaksHct116Tcf4V2.narrowPeak.hg18.nodups.gtf ~/genomic-data-misc/promoters/human_TSR_v1.Balwierz2009.gff.2000bp.gtf > ./wgEncodeYaleChIPseqPeaksHct116Tcf4V2.narrowPeak.NOTSS.hg18.gtf
	#
	#java GtfFilterNonoverlapping ./wgEncodeYaleChIPseqPeaksHct116Tcf4V2.narrowPeak.NOTSS.hg18.gtf  ~/genomic-data-misc/CGIs/Takai_Jones_plus_GG.merged.hg18.gtf > ./wgEncodeYaleChIPseqPeaksHct116Tcf4V2.narrowPeak.NOTSS.NOCGI.hg18.gtf
	#java GtfFilterNonoverlapping wgEncodeUwDnaseSeqPeaksRep1Hrcepic.narrowPeak.hg18.nodups.NOTSS.NOCGI.gtf wgEncodeUwDnaseSeqPeaksRep1Caco2.chr1.narrowPeak.hg18.nodups.NOTSS.NOCGI.gtf > wgEncodeUwDnaseSeqPeaksRep1.InHrecepic.NotInCaco2.NOTSS.NOCGI.gtf



	## Transform it to gtf.  If it's already GTF , skip this
	#my ($fh, $filename) = tempfile("buildMySqlTablesByChrom.pl.XXXXXX");
	my $inFn;
	my $outFn;
	my $cmd;
	if ($peakFn =~ /.g[tf]f$/i)
	{
		$outFn = $peakFn;
	}
	else
	{	
		$inFn = $peakFn;
		$outFn = "${name}.hg18.nodups.gtf";
		$cmd = "encodePeakFilesToGtf.pl ${inFn}";
		runCmd($cmd,0);
	}
	buildMethFile($outFn, $centeredSize);
	
	# Get Rid of TSS
	$inFn = $outFn;
	$outFn = "${name}.hg18.nodups.NOTSS2KB.gtf";
	my $cmd = "java -d64 -Xmx6000m GtfFilterNonoverlapping  ${inFn}  ~/genomic-data-misc/promoters/knownGene-tss.hg18.2000bp_centered.gtf > ${outFn}";
	runCmd($cmd,0);
	buildMethFile($outFn, $centeredSize);
	
	# Get rid of CGI
	$inFn = $outFn;
	$outFn = "${name}.hg18.nodups.NOTSS2KB.NOCGI.gtf";
	my $cmd = "java -d64 -Xmx6000m GtfFilterNonoverlapping  ${inFn}   ~/genomic-data-misc/CGIs/Takai_Jones_plus_GG.merged.hg18.gtf > ${outFn}";
	runCmd($cmd,0);
	buildMethFile($outFn, $centeredSize);

}

sub buildMethFile
{
	my ($fn, $centeredSize) = @_;
	
	my $centeredSec = ($centeredSize) ? "-centeredSize ${centeredSize}" : "";
	my $cmd = "java -d64 -Xmx6000m edu.usc.epigenome.scripts.MethylDbToFeatStatFiles ${centeredSec} methylCGsRich_normal010310_ methylCGsRich_tumor011010_ methylCGsRich_IMR90_ methylCGsRich_H1_ , ${fn}";
	runCmd($cmd,0);
}

sub runCmd
{
	my ($cmd, $testOnly) = @_;
	print STDERR $cmd."\n";
	print STDERR `$cmd`."\n" unless ($testOnly);
}
