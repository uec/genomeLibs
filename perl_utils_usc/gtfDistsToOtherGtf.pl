#!/usr/bin/env perl

use strict;
use List::Util qw[min max];
use File::Temp qw/tempfile/;

my $USAGE = "USAGE: gtfDistsToOtherGtf.pl from.gtf to.gtf > out.csv";


die "$USAGE\n" unless (scalar(@ARGV) == 2);

my ($fromFn, $toFn) = @ARGV;


# Run the java core
my ($tempfh, $tempfn) = tempfile();
close($tempfh);
my $cmd = "java NearestFeatures $fromFn $toFn > $tempfn";
print STDERR "Running: $cmd\n";
my $out = `$cmd`;

die "Can't read $tempfn\n" unless (open(F,$tempfn));


my $nlines = 0;
while (my $line = <F>)
{
	$nlines++;
	
	chomp $line;
	my @fldsOrig = split(/,/,$line);
	my ($dist,$fName,$fChr,$fSc,$fS,$fE,$fStr,$fType,$tName,$tChr,$tSc,$tS,$tE,$tStr,$tType) = @fldsOrig;
	
	# Determine more fine grained dist if it's inside and single point. (Make it negative)
	# The interpretation for non single points is more complicated
	my $singlePoint = ($fS==$fE);
	if ($singlePoint && ($dist==0))
	{
		my $distToLeft = $fS - $tS;
		my $distToRight = $tE - $fS;
		$dist = -1 * (min($distToLeft,$distToRight));
		#	print STDERR "L=$distToLeft\tR=$distToRight\t$dist\n";
	}
		
	#	print STDERR "$dist\t$fName ($fChr:$fS-$fE) $fSc -> $tName ($tChr:$tS-$tE) $tSc \n";

	# Reverse 0 and 1
	$fldsOrig[1] = $dist;
	$fldsOrig[0] = $fName;

	print join(",", qw/fName dist fChr fScore fStart fEnd fStrand fType tName tChr tScore tStart tEnd tStrand tType/)."\n" if ($nlines==1); # Header
	print join(",", @fldsOrig)."\n";
		
	
	
} 
close(F);
unlink($tempfn);



