#!/usr/bin/env perl

use File::Basename;
use File::Temp qw/ tempfile tempdir /;
use strict;
use Getopt::Long;
use File::Spec;

my $uecgatk = "$SOFTWAREROOT/uecgatk/2012-5-30/uecgatk.pl";
my $USAGE = "bamToElementEnrichmentInserts.pl [-distUpstream 1000] file.bam elements.bed output.txt";
my $TEMPPREFIX = "MATCHEDBEDINS";
my $TEMPDIR = ".";
my @REFS = ( "$GENOMEROOT/genomes/hg19_rCRSchrm/hg19_rCRSchrm.fa", 
	   "$GENOMEROOT/genomes/encode_hg19_mf/female.hg19.fa", 
	   "$GENOMEROOT/genomes/hg18_unmasked/hg18_unmasked.plusContam.fa");
@REFS = ( "$GENOMEROOT/genomes/hg19_rCRSchrm/hg19_rCRSchrm.fa");

my $bedhg18 = "$GENOMEROOT/CGIs/Takai_Jones_from_Fei_122007.fixed.hg18.PROMOTERONLY.oriented.bed";

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

	my ($fh, $outfn) = tempfile( "${TEMPPREFIX}counts.XXXXXX" , DIR => $TEMPDIR);
	close($fh);
	unlink($outfn);

	my $bambase = basename($bam,qw/\.bam/);
	my $bedbase = basename($bed,qw/\.bed/);
	my $tempbase = basename($outfn);
	$outfn = "${bambase}.${bedbase}.${tempbase}";


	my $cmd  = "$uecgatk -T ReadGroupProperties -R ${ref} -I ${bam} --intervals ${bed} -o ${outfn}";
	print STDERR "${cmd}...\n";
	`$cmd`;
	
	my $stdv = 0;
	my $medianTotal = 0;
	my $nMedians = 0;

	my $flds = readGatkReport($outfn);
	my $sizes = $flds->{"median.insert.size"} || [];

	foreach my $s (@$sizes)
	{
#	    print STDERR "size=$s\n";
	    $medianTotal += $s;
	    $nMedians++;
	}
	close (OUTF);

	my $mean = ($nMedians==0) ? 0 : ($medianTotal/$nMedians);
#	print STDERR "stdv=$stdv, mean=$mean ($bam)\n";
	return ($mean,$stdv);
}

sub readGatkReport
{
    my ($fn) = @_;

    die "Can't read ${fn}\n" unless (open(OUTF,$fn));
    my $colStarts = [];
    my $colEnds = [];
    my $fldNames = [];
    my $nCols = 0;
    my $out = {};
    while (my $line = <OUTF>)
    {
	chomp $line;
	if ($line =~ /^\#/)
	{
	}
	elsif ($line =~ /^\s*$/)
	{
	}
	elsif (!@$fldNames)
	{
	    my $lastE = -1;
	    while ($line =~ /(\S+\s*)/g) {
		my $fldName = $1;
		my $len = length($fldName);
		my $curS = $lastE + 1;
		my $curE = $curS + $len - 1;
		$fldName =~ s/\s*$//g;

		push(@$fldNames, $fldName);
		push(@$colStarts, $curS);
		push(@$colEnds, $curE);
		$out->{$fldName} = [];
		$nCols++;
		$lastE = $curE;
#		print STDERR "Found header col $nCols, ($fldName) \n";
	    }
	    # No header yet
	    /(\S+\s+)+/;
	}
	else
	{
	    # Data line. Extract fixed widths
	    for (my $i = 0; $i < $nCols; $i++)
	    {
		my $fldName = @${fldNames}[$i];
		my $fldList = $out->{$fldName};
		my $s = @${colStarts}[$i];
		my $e = @${colEnds}[$i];
		my $val = substr($line,$s,$e-$s+1);
		$val =~ s/\s*$//g;
		$val =~ s/^\s*//g;
#		print STDERR sprintf("Found data col (%s:%d-%d): \"%s\"\n", $fldName, $s,$e,$val);
		push(@$fldList,$val);
		$out->{$fldName} = $fldList;
	    }
	}
    }

    return $out;
}

sub makeMatchedElTempFile
{
	my ($inbed, $distUpstream) = @_;
	my $outbed = 0;
	my $fh = 0;
	
    print STDERR "CREATING temp file ${outbed} ($inbed, $distUpstream)\n";
    ($fh, $outbed) = tempfile( "${TEMPPREFIX}random.XXXXXX" , DIR => "/tmp");
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

