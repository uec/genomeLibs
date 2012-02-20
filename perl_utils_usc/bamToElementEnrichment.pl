#!/usr/bin/env perl

use File::Basename;
use File::Temp qw/ tempfile tempdir /;
use strict;
use Getopt::Long;

my $USAGE = "bamToElementEnrichment.pl [-distUpstream 1000] file.bam elements.bed ...";
my $TEMPPREFIX = "MATCHEDBED";
#my $SAMTOOLS = "/Users/benb/bin/samtools";
#my $COVERAGE = "/Users/benb/bin/BEDTools-Version-2.15.0/bin/bedtools coverage";
my $MEMPERJOB = "3995";
my $REF_FN_MAP = { "hg18" => "~/genomes/hg18_unmasked/hg18_unmasked.plusContam.fa",
	"hg19" => "~/genomes/hg19_rCRSchrm/hg19_rCRSchrm.fa" };
my $jarpath = "~/Java";
my $classpath = "${jarpath}/uecgatk/bin:${jarpath}/uecgatk/lib/GenomeAnalysisTK.jar:${jarpath}/uecgatk/lib/tuple.jar:${jarpath}/uecgatk/lib/commons-math-2.2.jar:" .
	"${jarpath}/uecgatk/lib/genomeLibs.jar:${jarpath}/uecgatk/lib/StingUtils.jar:${jarpath}/uecgatk/lib/UscKeck.jar";


my $distUpstream = 1000;
my $minMapq = 20;
GetOptions ('distUpstream=i', \$distUpstream, 'minq=i'=>\$minMapq) || die "$USAGE\n";

# Input params
die "$USAGE\n" unless (@ARGV==2);
my ($inbam, $inbed) = @ARGV;

# Get reference.
my $refFn = getRefFn($inbam);

# Run gatk counter on main file.
my ($amean, $amedian) = getCounts($inbam, $inbed, $minMapq);

# Run gatk counter on matched bed file
my ($matchedbed) = makeMatchedElTempFile($inbed, $distUpstream);
my ($bmean, $bmedian) = getCounts($inbam, $matchedbed, $minMapq);

my $ratio = $amean/$bmean;
print sprintf("RatioOfMeans=%0.3f\tbam=%s\tbed=%s\tdistUpstream=%d\n",$ratio,$inbam,$inbed,$distUpstream);

# Clean up
unlink($matchedbed);


# - - - - Functions

sub getRefFn
{
	my ($bam) = @_;
	my $genome = "hg19";
	if ($bam =~ /hg18/i)
	{
		$genome = "hg18";
	}
	elsif ($bam =~ /hg19/i)
	{
		$genome = "hg19";
	}
	else
	{
		print STDERR "Can't figure out reference genome from BAM filename: $bam\n";
		die;
	}
	
	return $REF_FN_MAP->{$genome};
}

sub getCounts
{
	my ($bam, $bed, $minMapq) = @_;
	my $cmd  = "";
	
	$cmd .= "export CLASSPATH=${classpath}";
	#$cmd .= "${SAMTOOLS} depth -Q ${minMapq} -b ${bed} ${bam}";
	my $cpus = 1;
	$cmd .= " ; java -Xmx${MEMPERJOB}m org.broadinstitute.sting.gatk.CommandLineGATK -T CoverageDepth -R ${refFn} -I ${bam} --intervals ${bed} -et NO_ET -nt ${cpus}"; # 
	
#	$cmd .= "${COVERAGE} -abam ${bam} -b ${bed} -hist | grep all";

	print STDERR "${cmd}...\n";
	my $out = `$cmd`;
	
	my $median = 0;
	my $mean = 0;
	LINE: foreach my $line (split(/\n/,$out))
	{
#		my ($chrom, $depth, $count, $total, @rest) = split(/\t/,$line);
#		$seenSoFar += $count;
#		my $fracSoFar = $seenSoFar / $total;
#		
#		print STDERR sprintf("depth=%d, count=%d, fracSoFar=%0.1f\n",$depth,$count, $fracSoFar*100);
#		$median = $depth;
#
#		last LINE if ($fracSoFar > 0.5);

		if ($line =~ /mean=(.*)$/)
		{ 
			$mean = $1;
		}
		elsif ($line =~ /50\.0 percentile=(.*)$/)
		{
			$median = $1;
		}
	}

	
	print STDERR "median=$median, mean=$mean ($bam)\n";
	return ($mean,$median);
}

sub makeMatchedElTempFile
{
	my ($inbed, $distUpstream) = @_;
	my $outbed = 0;
	my $fh = 0;
	
    print STDERR "CREATING temp file ${outbed} ($inbed, $distUpstream)\n";
    ($fh, $outbed) = tempfile( "${TEMPPREFIX}.XXXXXX" , DIR => "/tmp");
    die unless (open(BEDF,$inbed));
    
    my ($negStrand, $posStrand); 
    BEDLINE: foreach my $line (<BEDF>)
    {
    	chomp $line;
    	my @f = split(/\t/,$line);
    	my $size = $f[2] - $f[1];
    	
    	if ($f[5] eq '-')
    	{
    		$f[1] = $f[2] + $distUpstream;
    		$f[2] = $f[1] + $size;
    		$negStrand++;
    	}
    	else
    	{
    		$f[2] = $f[1] - $distUpstream;
    		$f[1] = $f[2] - $size;
    		$posStrand++;
    	}
    	print $fh join("\t",@f)."\n";
    	
    }
   	print STDERR sprintf("Saw %d pos strand and %d neg strand\n",$posStrand,$negStrand);
	close($fh);
	
	# Stupid GATK won't accept it unless it ends with bed
	my $newoutbed = "${outbed}.bed";
	`mv $outbed $newoutbed`;
	
	return $newoutbed;
}


