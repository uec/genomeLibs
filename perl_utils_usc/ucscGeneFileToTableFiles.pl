#!/usr/bin/perl

use strict;
use File::Basename;
use Getopt::Long;

my $USAGE = "ucscGeneFileToTableFiles.pl --delim \",\" (default tabl) --outPrefix features --desc knownGenes < input.txt";

my $delim = "\t";
my $outPrefix = "features";
my $desc = "";
GetOptions ('desc=s' => \$desc, 'outPrefix=s' => \$outPrefix, 'delimx=s' => \$delim) || die "$USAGE\n";

die "$USAGE\n" unless (@ARGV==0);


# Setup output files
$::openFhs = {};

my $fieldmap = {};
my $lineNum = 0;
my $seen = {};
LINE: while (my $line = <STDIN>)
{
	$lineNum++;
	print STDERR "On line $lineNum\n" if (($lineNum%1000)==0);
	
	chomp $line;
	my @f = split($delim, $line);
	
	# Check for valid line
	if ($line =~ /knownGene\.name/)
	{
		print STDERR "Found header\n";
		# Header
		for (my $i = 0; $i < scalar(@f); $i++)
		{
			$fieldmap->{$f[$i]} = $i;
		}
		next LINE;
	}
	else
	{
		die "No header line!\n" unless (scalar(keys(%$fieldmap))>0); 
	}
	
	next LINE if ($line =~ /^\s+$/);
	#print STDERR "Line $lineNum num flds: " . scalar(@f) . "\n";
	die "Illegal line: $line\n" unless ((@f==15) || (@f==14));
		
	# Parse line
	my $chr = $f[$fieldmap->{"hg18.knownGene.chrom"}];
	next LINE if ($chr =~ /random/);
	next LINE if ($chr =~ /hap/);

	my $strand = $f[$fieldmap->{"hg18.knownGene.strand"}];
	my $txS = $f[$fieldmap->{"hg18.knownGene.txStart"}];
	my $txE = $f[$fieldmap->{"hg18.knownGene.txEnd"}];
	my $exonStartsStr = $f[$fieldmap->{"hg18.knownGene.exonStarts"}];
	my $exonEndsStr = $f[$fieldmap->{"hg18.knownGene.exonEnds"}];
	my $geneSymbol = $f[$fieldmap->{"hg18.kgXref.geneSymbol"}];
	my $refseqId = $f[$fieldmap->{"hg18.kgXref.refseq"}];
	
	
	
	my $refseqVers = 0;
	if ($refseqId =~ /^(.*)(\.[0-9]+)$/)
	{
		$refseqId = $1;
		$refseqVers = $2;
	}
	
	# Start a template
	my @outFlds = (0, 0, 0, $strand, "0.0", $refseqId, $refseqVers, $geneSymbol);
	
	# simple ones
	writeRow($outPrefix, $desc, $chr, "tx", $txS, $txE, $strand, "0.0", $refseqId, $refseqVers, $geneSymbol);
	my $tss = ($strand eq '-') ? $txE : $txS;
	writeRow($outPrefix, $desc, $chr, "tss", $tss, $tss, $strand, "0.0", $refseqId, $refseqVers, $geneSymbol);
	
	# Do some extra ones for TSS
	writeRow($outPrefix, $desc, $chr, "tss_2kb_flank", $tss-2000, $tss+2000, $strand, "0.0", $refseqId, $refseqVers, $geneSymbol);
	writeRow($outPrefix, $desc, $chr, "tss_1kb_flank", $tss-1000, $tss+1000, $strand, "0.0", $refseqId, $refseqVers, $geneSymbol);
	writeRow($outPrefix, $desc, $chr, "tss_500bp_flank", $tss-500, $tss+500, $strand, "0.0", $refseqId, $refseqVers, $geneSymbol);
	

	# And exons
	my @exonStarts = split(/,/,$exonStartsStr);
	my @exonEnds = split(/,/,$exonEndsStr);
	die "Why do we have a different number of exon start than ends?\n${exonStartsStr}\n${exonEndsStr}\n"
		unless (scalar(@exonStarts) == scalar(@exonEnds));
	for (my $i = 0; $i < @exonStarts; $i++)
	{
		my $exonS = $exonStarts[$i];
		my $exonE = $exonEnds[$i];
		writeRow($outPrefix, $desc, $chr, "exon", $exonS, $exonE, $strand, "0.0", $refseqId, $refseqVers, $geneSymbol);
	}

}

# close files
foreach my $fh (values(%$::openFhs))
{
	close($fh);
}

# - - - - Functions

sub writeRow
{
	my ($prefix, $desc, $chr, $feat, $s, $e, $strand, $score, $refseqId, $refseqVers, $name) = @_;
	
	my $descSec = ($desc) ? ".${desc}" : "";
	my $fn = "${prefix}_${chr}${descSec}.txt";
	my $fh = $::openFhs->{$fn};
	if (!$fh)
	{
		die "Can't write to $fn\n" unless (open($fh,">$fn"));
		$::openFhs->{$fn} = $fh; 
	}
	
	print $fh join("\t",$feat, $s, $e, $strand, $score, 
		($refseqId) ? $refseqId : "NULL", 
		($refseqVers) ? $refseqVers : "0", $name);
	print $fh "\n";
}