#!/usr/bin/env perl

use File::Basename;
use File::Temp qw/ tempfile tempdir /;
use strict;
use Getopt::Long;

my $USAGE = "bamToElementEnrichment.pl [-distUpstream 1000] file.bam elements.bed ...";
my $TEMPPREFIX = "MATCHEDBED";
#my $SAMTOOLS = "/Users/benb/bin/samtools";
my $COVERAGE = "/Users/benb/bin/BEDTools-Version-2.15.0/bin/bedtools coverage";
my $MEMPERJOB = "3995";

#my $jarpath = "~/Java";
#my $classpath = "${jarpath}/uecgatk/bin:${jarpath}/uecgatk/lib/GenomeAnalysisTK.jar:${jarpath}/uecgatk/lib/tuple.jar:${jarpath}/uecgatk/bin:${jarpath}/uecgatk/lib/GenomeAnalysisTK.jar:" .
#	"${jarpath}/uecgatk/lib/genomeLibs.jar:${jarpath}/uecgatk/lib/StingUtils.jar:${jarpath}/uecgatk/lib/UscKeck.jar";
#my $refFn = "~/uec-00/genomes/hg18_unmasked/hg18_unmasked.plusContam.fa";

my $distUpstream = 1000;
my $minMapq = 20;
GetOptions ('distUpstream=i', \$distUpstream, 'minq=i'=>\$minMapq) || die "$USAGE\n";

# Input params
die "$USAGE\n" unless (@ARGV==2);
my ($inbam, $inbed) = @ARGV;


# Run gatk counter on main file.
my ($a) = getCounts($inbam, $inbed, $minMapq);

# Run gatk counter on matched bed file
my ($matchedbed) = makeMatchedElTempFile($inbed, $distUpstream);
my ($b) = getCounts($inbam, $matchedbed, $minMapq);


# Clean up
#unlink($matchedbed);


sub getCounts
{
	my ($bam, $bed, $minMapq) = @_;
	my $cmd  = "";
#	$cmd .= "export CLASSPATH=${classpath}";
#	#$cmd .= "${SAMTOOLS} depth -Q ${minMapq} -b ${bed} ${bam}";
#	my $cpus = 1;
#	$cmd .= " ; java -Xmx${MEMPERJOB}m org.broadinstitute.sting.gatk.CommandLineGATK -T DepthOfCoverage -R ${refFn} -I ${bam} -mmq ${minMapq}  --intervals ${bed} -et NO_ET -nt ${cpus}";
	
	$cmd .= "${COVERAGE} -abam ${bam} -b ${bed} -hist | grep all";
	print STDERR "${cmd}...\n";
	my $out = `$cmd`;
	
	my $seenSoFar = 0;
	my $median = 0;
	LINE: foreach my $line (split(/\n/,$out))
	{
		my ($chrom, $depth, $count, $total, @rest) = split(/\t/,$line);
		$seenSoFar += $count;
		my $fracSoFar = $seenSoFar / $total;
		
		print STDERR sprintf("depth=%d, count=%d, fracSoFar=%0.1f\n",$depth,$count, $fracSoFar*100);
		$median = $depth;

		last LINE if ($fracSoFar > 0.5);
	}
	
	print STDERR "median=$median\n";
	return $median;
}

sub makeMatchedElTempFile
{
	my ($inbed, $distUpstream) = @_;
	my $outbed = 0;
	my $fh = 0;
	
    print STDERR "CREATING temp file ${outbed} ($inbed, $distUpstream)\n";
    ($fh, $outbed) = tempfile( "${TEMPPREFIX}.XXXXXX" , DIR => "/tmp");
    die unless (open(BEDF,$inbed));
    
    BEDLINE: foreach my $line (<BEDF>)
    {
    	chomp $line;
    	my @f = split(/\t/,$line);
    	my $size = $f[2] - $f[1];
    	
    	if ($f[5] eq '-')
    	{
    		$f[1] = $f[2] + $distUpstream;
    		$f[2] = $f[1] + $size;
    	}
    	else
    	{
    		$f[2] = $f[1] - $distUpstream;
    		$f[1] = $f[2] - $size;
    	}
    	print $fh join("\t",@f);
    }


	close($fh);
	return $outbed;
}


