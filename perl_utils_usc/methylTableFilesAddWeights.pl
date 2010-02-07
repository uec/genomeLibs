#!/usr/bin/env perl

use strict;
use File::Basename;
use Getopt::Long;
use File::Temp qw/ tempfile /;

my $USAGE = "methylTableFilesAddWeights.pl -cpgLocFn cpgs.cse methylCpgsRich_IMR90_chr1.tab file2.tab file3.tab ...";

my $cpgLocFn = "~/genomic-data-misc/CpgsAll.hg18.cse";
GetOptions ('cpgLocFn=s'=>\$cpgLocFn)|| die "$USAGE\n";

die "$USAGE\n" unless (@ARGV>0);

# Read in weights from cpg file
print STDERR "Reading weights from $cpgLocFn\n";
my $weights = readWeightsFromCse($cpgLocFn);
print STDERR "Found " . scalar(keys(%$weights)) . " Cpg weights\n";

FILE: foreach my $origFn (@ARGV)
{
    my ($name, $path, $suf) = fileparse($origFn, qr/\.[^.]*/); # qr/\.txt/); #

	if ($name !~ /chr([0-9xy]+)/i)
	{
		print STDERR "File $name does not contain \"chr\*\" .. skipping\n";
		next FILE;
	}
	my $chrNum = $1;
	$chrNum = 23 if (uc($chrNum) eq "X");
	$chrNum = 24 if (uc($chrNum) eq "Y");
	$chrNum = 25 if (uc($chrNum) eq "M");

	# We want the output file to have the same name as the input (because
	# this matches the database name).  So we move the original to a temp.
	my ( $newfh, $newFn ) = tempfile( "${name}XXXXXX", DIR => "/tmp" );
	close($newfh);
	`cp $origFn $newFn`;  # Change to mv	

	# Write to original file
	die "Can't write to $origFn\n" unless (open(W, ">$origFn"));
	
	# Read from temp file
	die "Can't read from $newFn\n" unless (open(R, $newFn));

	LINE: while (my $line = <R>)
	{
		chomp $line;
		my @f = split(/\t/,$line);
		
		my $pos = $f[0];
		$pos-- if ($f[1] eq "-"); # Pos should always be the forward strand C
		my $key = "${chrNum}.${pos}";
		my $weight = $weights->{$key};
		
		if (!$weight)
		{
			print STDERR  "Can't find weight for CpG at $key\n";
			$weight = 1;
		} 
		
		my $replace = 0;
		$replace = 1 if ($f[11] =~ /^[0-9]+$/i);
		
		splice(@f, 11, ($replace) ? 1 : 0, $weight);
		print W join("\t",@f)."\n";
	}
	
	close(R);
	close(W); 
}

sub runCmd
{
	my ($cmd, $testOnly) = @_;
	print STDERR $cmd."\n";
	print STDERR `$cmd`."\n" unless ($testOnly);
}

sub readWeightsFromCse
{
	my ($fn) = @_;
	
	die "Can't read $fn\n" unless (open(F, $fn));
	my $weights = {};
	my $lastPos = -1;
	my $lastLastPos = -1;
	my $lastChr = -1;
	my $lastLastChr = -1;
	my $cpgNum = 0;
	while (my $line = <F>)
	{
		$cpgNum++;
		if (($cpgNum % 1000000)==0)
		{
			print STDERR "On Cpg $cpgNum (chr=$lastChr, pos=$lastPos)\n";
		}
		
		chomp $line;
		my @f = split(/,/,$line);
		my ($chr, $pos, $nextpos) = @f;
		
		# Add 1 to make it 1-based
		$pos++; $nextpos++;

		# Setup up stock cases, which are not correct at the beginning
		# and ends of chromosomes.
		my $windS = int(($lastLastPos + $lastPos)/2);
		my $windE = int(($lastPos + $pos)/2);
		
		# Special cases for the last position on one chrom and 
		# the first of the next one
		my $include = 0;
		if (($lastPos<0) && ($lastLastPos<0))
	 	{
	 		# On the first one
	 		#print STDERR "On first one\n";
	 		$include = 0;
	 	}
		elsif ($lastChr ne $chr)
		{
			# "last" is the last one on prev chrom
			$windE = $lastPos;
			$include = 1;
		}
		elsif ($lastLastChr ne $lastChr)
		{
			# "last" is the first one on new chrom
			$windS = $lastPos;
			$include = 1;
		}
# Case 2 already takes care of this	 	
#	 	elsif (($lastPos>=0) && ($lastLastPos<0))
#	 	{
#	 		# On the second to last one
#	 		$windS = $lastPos;
#	 		$include = 1;
#	 	}
	 	else
	 	{
	 		# We're in the middle of chrom
	 		#print STDERR "In normal position\n";
	 		$include = 1;
	 	}
	 	
	 	# Add to weights
	 	if ($include)
	 	{
	 		my $weight = $windE-$windS+1;
	 		my $key = "${lastChr}.${lastPos}";
	 		#print STDERR "\tweights->{$key}=$weight (${windS}-${windE})\n";
	 		$weights->{$key} = $weight;
	 	}
		
		# Increment
		$lastLastPos = $lastPos;
		$lastPos = $pos;
		$lastLastChr = $lastChr;
		$lastChr = $chr;
	}
	close(F);
	
	# We still have to do the last one in the genome (now lastPos and lastChr)
	my $windS = int(($lastLastPos + $lastPos)/2);
	my $windE = $lastPos;
	my $weight = $windE-$windS+1;
	my $key = "${lastChr}.${lastPos}";
	$weights->{$key} = $weight;
		
	
	return $weights;
}



