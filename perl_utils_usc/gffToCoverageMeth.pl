#!/usr/bin/env perl

use File::Basename;
use File::Temp qw/ tempfile tempdir /;
use strict;
use Getopt::Long;

my $USAGE = "gffToCoverageMeth.pl -bam a.bam file1.gff file2.gff -ref hg19_rCRSchrm.fa";

my $bam = "";
my $reffa = "~/genomes/hg19_rCRSchrm/hg19_rCRSchrm.fa";
GetOptions ('bam=s',\$bam,
	    'ref=s',\$reffa,
	    ) || die "$USAGE\n";
# 'compareCutoff=f',\$compareCutoff,'testCutoff=f',\$testCutoff, 
#print STDERR "doBare=${doBare}\tintervalFile=${intervalFile}\n";

# Input params
die "$USAGE\n" unless (@ARGV>=1);
my (@gffs) = @ARGV;

 my $bambase = basename($bam, qw/.bam/);

foreach my $gff (@gffs)
{
    my $base = basename($gff, qw/.g[tf]f/);
    my $gffout = "cvg.${base}.${bambase}.csv";

    print STDERR "$gff -> $gffout\n";
    unlink($gffout);

    die "Can't read $gff\n" unless (open(FIN,$gff));

    my $onLine = 0;
    while (my $line = <FIN>)
    {
	next if ($line =~ /^\s*track/);
	next if ($line =~ /^\s*\#/);
	chomp $line;

	my ($chr, $a, $b, $start, $end, $c, $d, $e, $f) = split(/\t/,$line);
	outputRegion($chr, $start, $end, $bam, $gffout,$onLine==0);
	$onLine++;
    }

    close(FIN);
}

sub outputRegion
{
    my ($chr, $s, $e, $bam, $outfn, $usehead) = @_;

    my $chrNum = $chr;
    $chrNum =~ s/X/23/ig;
    $chrNum =~ s/Y/24/ig;
    $chrNum =~ s/M/25/ig;
    $chrNum =~ s/chr//g;

    my $loc = sprintf("%s:%d-%d",$chr,$s,$e);

    # Make a temp file
#    my  ($tempfh, $tempfn) = tempfile(DIR => "./");
   my  ($tempfh, $tempfn) = tempfile();
    close($tempfh);
    my $base = basename($tempfn);

    my $headsec = ($usehead) ? "--header" : "";
 
    my $cmd = "java -Xmx5495m org.broadinstitute.sting.gatk.CommandLineGATK -T MethLevelAveragesCoverage -R ${reffa} -I ${bam} --iupacPatterns NCG  -nt 1  -L ${loc} --out ${tempfn} ${headsec}";
    print STDERR "Running $cmd\n";
    `$cmd`;

    die "Can't open ${tempfn}\n" unless (open(TMP,$tempfn));
    die "Can't write to ${outfn}\n" unless (open(OUTF,">>${outfn}"));
    while (my $line = <TMP>)
    {
	chomp $line;
	if ($line =~ /NumReads/)
	{
	    # Header
	    $line = "chr,start,end,${line}\n";
	}
	else
	{
	    $line = "${chrNum},${s},${e},${line}\n";
	}
	print OUTF $line;
    }
    close(OUTF);
    close(TMP);
    unlink($tempfn);
}
