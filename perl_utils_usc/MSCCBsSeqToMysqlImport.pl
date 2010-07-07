#!/usr/bin/env perl

use strict;
use File::Basename;
use Getopt::Long;
use POSIX;

my $USAGE = "MSCCBsSeqToMysqlImport.pl --outPrefix mnase < input.txt";

my $delim = "\t";
my $outPrefix = "features";
my $desc = "";
GetOptions ('outPrefix=s' => \$outPrefix) || die "$USAGE\n";

die "$USAGE\n" unless (@ARGV==0);


# Setup output files
$::openFhs = {};

#my $fieldmap = {};
my $lineNum = 0;
#my $seen = {};
#my $numDups = 0;
#my $lastKey = "";
#my $lastChr = 0;
my @countsTotal = (0,0,0);
LINE: while (my $line = <STDIN>)
{
	$lineNum++;
	print STDERR "On line $lineNum\n" if (($lineNum%100000)==0);
	
	chomp $line;
	my @f = split(/\t/, $line);
	
	# Check for valid line
	if ((@f != 8) || ($line !~ /^chr/))
	{
		print STDERR "Unrecognized Cokus format! $line\n";
		next LINE;
	} 
	
	my ($chr, $pos, $strand) = @f[0..2];
	next LINE if ($chr =~ /_/);  # Random and haplotypes
	
	# They put them as the same one.  I double checked this matches my coordinate system.
	$pos++ if ($strand eq '-');

	my @counts = @f[3..5];
	foreach my $i (0..2) {$countsTotal[$i] += $counts[$i]} 
	my $count = $counts[0] + $counts[1] + $counts[2];
	my $meth = countToMeth($count);
	
	#printf STDOUT ("%s\t\t\t%d\t%0.2f\n", $line, $count, $meth) ;
	

	writeRow($outPrefix, $chr, $pos, $strand, 100, sprintf("%0.0f",$meth*100));

}



print STDERR "Counts total: ".join(", ", @countsTotal)."\n";

# close files
foreach my $fh (values(%$::openFhs))
{
	close($fh);
}

# - - - - Functions

sub countToMeth
{
	my ($count) = @_;
	my $meth = (-0.1124 * $count) + 1;
	$meth = 0 if ($meth<0);
	return $meth;
}

sub writeRow
{
	my ($prefix, $chr, $pos, $strand, $total, $totalMeth) = @_;
	my $fn = "${prefix}_${chr}.txt";
	my $fh = $::openFhs->{$fn};
	if (!$fh)
	{
		die "Can't write to $fn\n" unless (open($fh,">$fn"));
		$::openFhs->{$fn} = $fh; 
	}
	
	print $fh join("\t",$pos, $strand, $total, $totalMeth, 0, ($total-$totalMeth), 0, 0, 0, 0, 0, 0);
	print $fh "\n";
}

