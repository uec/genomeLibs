#!/usr/bin/env perl

use File::Basename;
use File::Temp qw/ tempfile tempdir /;
use strict;
use Getopt::Long;
use File::Spec;

my $PICARD = "/home/uec-00/shared/production/software/picard/default/";
my $JAVA = "/home/uec-00/shared/production/software/java/default/bin/java";
my $uecgatk = "/home/uec-00/shared/production/software/uecgatk/default/uecgatk.pl";
my $USAGE = "bamToElementEnrichment.pl [-distUpstream 1000] file.bam elements.bed output.txt";
my $TEMPPREFIX = "MATCHEDBED";
my @REFS = ( "/home/uec-00/shared/production/genomes/hg19_rCRSchrm/hg19_rCRSchrm.fa", 
	   "/home/uec-00/shared/production/genomes/encode_hg19_mf/female.hg19.fa", 
	   "/home/uec-00/shared/production//genomes/hg18_unmasked/hg18_unmasked.plusContam.fa",
	   "/home/uec-00/shared/production//genomes/mm9_unmasked/mm9_unmasked.a",
	   "/home/uec-00/shared/production/genomes/mm10/mm10.fa");

my $bedhg18 = "/home/rcf-40/bberman/tumor/genomic-data-misc/CGIs/Takai_Jones_from_Fei_122007.fixed.hg18.PROMOTERONLY.oriented.bed";
my $mm9bed = "/home/uec-00/shared/production/genomic-data-misc/CpG_islands/mm9cpgisland.bed";
my $mm10bed = "/home/uec-00/shared/production/genomic-data-misc/CpG_islands/mm10cpgisland.bed";

my $distUpstream = 1000;
my $minMapq = 20;
GetOptions ('distUpstream=i', \$distUpstream, 'minq=i'=>\$minMapq) || die "$USAGE\n";

# Input params
die "$USAGE\n" unless (@ARGV==3);
my ($inbam, $inbed,$outfile) = @ARGV;


if (calculateRatio($inbam, $inbed, $minMapq) == 0)
{
	my $basebam = basename($inbam);
	my $workdir = $basebam . time();

	#mkdir($workdir);
	#chdir($workdir);
	#system("$JAVA -Xmx10g -jar $PICARD/SortSam.jar CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT INPUT=$inbam OUTPUT=$basebam\.coordsort.bam SORT_ORDER=coordinate");
	#system("cp $basebam\.coordsort.bai $basebam\.coordsort.bam.bai");
	#calculateRatio("$basebam\.coordsort.bam", $inbed, $minMapq);
}

exit;


sub calculateRatio
{
	my ($bam, $bed, $minMapq) = @_;
	for my $ref (@REFS)
	{

		$bed = $bedhg18 if($ref =~ /hg18/);
		$bed = $mm10bed if($ref =~ /mm10/);
		$bed = $mm9bed if($ref =~ /mm9/);
		# Run gatk counter on main file.
		my ($amean, $astdv) = getCounts($bam, $bed, $minMapq, $ref);
		if($amean)
		{

			# Run gatk counter on matched bed file
			my ($matchedbed) = makeMatchedElTempFile($bed, $distUpstream);
			my ($bmean, $bstdv) = getCounts($bam, $matchedbed, $minMapq, $ref);
			my $ratio = $amean/$bmean;
			open(OUT,"> $outfile");
			print OUT sprintf("RatioOfMeans=%0.3f\tamean=%0.3f\tbmean=%0.3f\tastdv=%0.3f\tbstdv=%0.3f\tbam=%s\tbed=%s\tref=%s\tdistUpstream=%d\n",$ratio,$amean, $bmean, $astdv, $bstdv, File::Spec->rel2abs($inbam),$bed,$ref,$distUpstream);
			close OUT;

			# Clean up
			unlink($matchedbed);
			return 1;
		}
	}
	return 0;
}


# - - - - Functions

sub getCounts
{
	my ($bam, $bed, $minMapq,$ref) = @_;
	my $cmd  = "$uecgatk -T CoverageDepth -R ${ref} -I ${bam} --intervals ${bed} ";
	print STDERR "${cmd}...\n";
	my $out = `$cmd`;
	
	my $stdv = 0;
	my $mean = 0;
	foreach my $line (split(/\n/,$out))
	{
		if ($line =~ /mean=(.*)$/)
		{ 
			$mean = $1;
		}
		elsif ($line =~ /std dev=(.*)$/)
		{
			$stdv = $1;
		}
	}

	
	print STDERR "stdv=$stdv, mean=$mean ($bam)\n";
	return ($mean,$stdv);
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
    	print $fh join("\t",@f)."\n" if ($f[1] > 0 && $f[2] > 0);
    	
    }
   	print STDERR sprintf("Saw %d pos strand and %d neg strand\n",$posStrand,$negStrand);
	close($fh);
	
	# Stupid GATK won't accept it unless it ends with bed
	my $newoutbed = "${outbed}.bed";
	`mv $outbed $newoutbed`;
	
	return $newoutbed;
}

