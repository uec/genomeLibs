#!/usr/bin/env perl

use strict;
use File::Basename;
use Getopt::Long;
use File::Temp qw/ tempfile /;

my $USAGE = "FeatTableIntersection.pl [--db cr] [--testOnly] outputFeatType filterFeatType";
my @CHROMS = map {"chr" . $_} (1..22,'X','Y','M');


my $db = "cr";
my $testOnly = 0;
GetOptions ('db=s'=>\$db, 'testOnly!',\$testOnly) || die "$USAGE\n";
die "$USAGE\n" unless (@ARGV==2);
my ($outputFeatType, $filtFeatType) = @ARGV;

my $foundHeader = 0;
foreach my $chr (@CHROMS)
{
	my $groupBy = "GROUP BY chromPosStart,chromPosEnd";
#	$groupBy = "GROUP BY prim.name";
	my $sql = "select prim.* from features_${chr} prim, features_${chr} filt where " .
	 	"NOT ((prim.chromPosStart>filt.chromPosEnd) OR (prim.chromPosEnd<filt.chromPosStart)) AND " . 
	 	"(prim.featType='${outputFeatType}') AND (filt.featType='${filtFeatType}') " .
	 	"${groupBy} ORDER BY prim.chromPosStart;";
	#print STDERR "Chrom $sql\n";
	
	my $cmd = "echo \"$sql\" | mysql ${db}";
	my $output = runCmd($cmd, $testOnly);
	
	LINE: foreach my $line (split(/\n/,$output))
	{
		#my ($featType, $s, $e, $strand, $score, $refSeqId, $refSeqVers, $symbol) = split(/\t/,$line);
		if ($line =~ /featType/) # Header
		{
			print "chrom\t$line\n" if (!$foundHeader);
			$foundHeader = 1; 
		}
		else
		{
			print "$chr\t$line\n";
		}
	}
	
}


sub runCmd
{
	my ($cmd, $testOnly) = @_;
	print STDERR $cmd."\n";
	my $out = 0;
	if (!$testOnly)
	{
		$out = `$cmd`;
	}
	return $out;
}
