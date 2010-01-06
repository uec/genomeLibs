#!/usr/bin/perl

use strict;
use File::Basename;
use Getopt::Long;

my $USAGE = "methylDbIntersectFeatsWithExpression.pl outPrefix";

my $DB = "cr";
my @chroms = ("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", 
	       "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", 
	       "chr19", "chr20", "chr21", "chr22", "chrX", "chrY", "chrM");
#@chroms = ("chr11");

# Rank fields go from lowest expression to highest
my $filters = 
{
	tumHigh => "e.pctileTumor > 0.667",
	tumLow => "e.pctileTumor < 0.333",
	normHigh => "e.pctileNormal > 0.667",
	normLow => "e.pctileTumor < 0.333",
	tumUp => "e.ttest < 0.05 AND e.meanTumor > e.meanNormal",
	tumDown => "e.ttest < 0.05 AND e.meanTumor < e.meanNormal"
};


#GetOptions ('chromCol=i' => \$chromCol, 'omitChromFldInOutput' => \$omitChrFld, 'outPrefix=s' => \$outPrefix, 'delimx=s' => \$delim) || die "$USAGE\n";


die "$USAGE\n" unless (@ARGV>0);
my ($outprefix) = @ARGV;
$outprefix .= "_" unless ($outprefix =~ /_$/);
 
 
 foreach my $chr (@chroms)
 {
 	foreach my $filt (keys(%$filters))
 	{
  		my $whereSecFilt = $filters->{$filt};
  		my $whereSecJoin = "e.refseqId = f.refseqId";

		my $sql = "select concat(f.featType,\"_${filt}\"),f.chromPosStart, f.chromPosEnd, f.strand, f.score, f.refseqId, f.refseqVers, f.name ".
			" FROM infiniumExpr_${chr} e, features_${chr} f WHERE ${whereSecFilt} AND ${whereSecJoin};";
		print STDERR "$sql\n";

		my $outfn = "${outprefix}${chr}.${filt}.txt";
		my $cmd = "echo '${sql}' | mysql $DB > $outfn";
		print STDERR "$cmd\n";
		print STDERR `$cmd\n`;
		

	}
	
	my $cmd = "cat ${outprefix}*${chr}\.*txt > ${outprefix}${chr}.txt";
	print STDERR "$cmd\n";
	print STDERR `$cmd\n`;
 }