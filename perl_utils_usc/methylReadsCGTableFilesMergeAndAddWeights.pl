#!/usr/bin/env perl

use strict;
use File::Basename;
use Getopt::Long;
use File::Temp qw/ tempfile /;

my $NOWEIGHTS = 0;
my $USAGE = "methylReadsCGTableFilesMergeAndAddWeights.pl -cpgLocFn cpgs.cse -outPrefix methylReadCG methylReadCG_PCR1 methylReadCG_PCR2 ...";

my $cpgLocFn = "/Users/benb/genomic-data-misc/CpG_islands/CpgsAll.hg18.cse";
my $outPrefix = "methylReadCG";
GetOptions ('cpgLocFn=s'=>\$cpgLocFn, 'outPrefix=s'=>\$outPrefix) || die "$USAGE\n";

die "$USAGE\n" unless (@ARGV>0);

# Read in weights from cpg file
print STDERR "Reading weights from $cpgLocFn\n";
my $weights = {};
$weights = readWeightsFromCse($cpgLocFn) unless ($NOWEIGHTS);
print STDERR "Found " . scalar(keys(%$weights)) . " Cpg weights\n";


# Go through chroms
my @chrs = 1..25;
@chrs = (1);
my @prefixes = @ARGV;
foreach my $chrNum (@chrs)
{

	my $chrRowsH = {};

	# First pass just to check filenames
	my $foundFiles = fillChromHash($chrRowsH, $weights, $chrNum, \@prefixes, 1, 0);
	if (!$foundFiles)
	{
		print STDERR "Could not find one or more input files, quitting.\n";
		exit(1);
	}

	# Second pass to get chrom hash
	print STDERR "Starting second pass...\n";
	fillChromHash($chrRowsH, $weights, $chrNum, \@prefixes, 0, 0);
	print STDERR sprintf("Found %d entries\n",scalar(keys(%$chrRowsH)));

	# Now dump the actual rows
	print STDERR "Starting third pass...\n";
	fillChromHash($chrRowsH, $weights, $chrNum, \@prefixes, 0, 1);
			
#	# Then go through and output
#	outputHash($chrRowsH, $outPrefix, $chrNum);
}


#sub outputHash
#{
#	my ($chrRowsH, $outPrefix, $chrNum) = @_;
#	
#	my $chr = chrNumToChr($chrNum);
#
#	my $outfn = sprintf("%s_%s.txt", $outPrefix, $chr);
#	die "Can't write to $outPrefix\n" unless (open(W,">$outfn"));
#	foreach my $pos (sort(keys(%$chrRowsH)))
#	{
#		my $row = $chrRowsH->{$pos};
#		print W join("\t", @$row);
#		print W "\n";
#	}
#	
#	close(W);
#}

sub chrNumToChr
{
	my ($chrNum) = @_;
	
	my $chr = "chr${chrNum}";
	$chr = "chrX" if ($chrNum==23);
	$chr = "chrY" if ($chrNum==24);
	$chr = "chrM" if ($chrNum==25);
	return $chr;
}

# Returns true if it was able to read all the files
sub fillChromHash
{
	my ($chrRowsH, $weights, $chrNum, $prefixes, $fileCheckOnly, $dumpRows) = @_;

	my $chr = chrNumToChr($chrNum);
	
	my $outfn = sprintf("%s_%s.txt", $outPrefix, $chr);
	if ($dumpRows)
	{
		die "Can't write to $outPrefix\n" unless (open(W,">$outfn"));
	}
	
	
	my $out = 1;
	my $fileNum = 0;
	FILE: foreach my $prefix (@$prefixes)
	{
		
		# Make sure the file exists
		my $fn = sprintf("%s_%s.txt",$prefix,$chr);
		print STDERR "Looking for file ${fn}\n"; 
		
		my $foundFn = (-f $fn);
		$out &&= $foundFn;
		print STDERR "Can't find file $fn\n" unless ($foundFn);
					
		# We're done if it's just a file check
		next FILE if ($fileCheckOnly || !$foundFn);		
			


	

		# Read from temp file
		die "Can't read from $fn\n" unless (open(R, $fn));
	
		my $keysSeenThisFile = {};
	
		LINE: while (my $line = <R>)
		{
			chomp $line;
			my @f = split(/\t/,$line);
			
			
			my $pos = $f[0];
			#print STDERR "orig pos = $pos\n";
			my $weightPos = ($f[1] eq "-") ? ($pos-1) : $pos; # Pos should always be the forward strand C
			#print STDERR "first strand correction = $pos\n";
			my $weightKey = "${chrNum}.${weightPos}";
			my $weight = $weights->{$weightKey};
			
			# The H1 and IMR90 dumps didn't have negative strands, which we can detect
			# by lack of a weight entry
			if (!$NOWEIGHTS && !$weight)
			{
				if ($f[1] eq '+')
				{

					$weightPos--;
					#print STDERR "second strand correction = $pos\n";
					$weightKey = "${chrNum}.${weightPos}";
					$weight = $weights->{$weightKey};
					$f[1] = "-" if ($weight); # This is actually a negative strand one
				}
				
				if (!$weight)
				{
					print STDERR  "Can't find weight for CpG at $weightKey\n";
					$weight = 1;
				}
			} 
			
			# Make sure readIds don't conflict
			$f[2] += ($fileNum * 1E9);


			my $replace = 0;
			$replace = 1 if (($f[11] =~ /^[0-9]+$/i) || ($f[11] =~ /^nan$/i));
			
			print STDERR "replace=$replace, numFlds = " . scalar(@f) . "\n";
			if (scalar(@f)>=12)
			{
				splice(@f, 11, ($replace) ? 1 : 0, $weight);
			}
			else
			{
				for (my $i = scalar(@f); $i < 11; $i++)
				{
					$f[$i] = 0;
				}
				$f[11] = $weight;
			}
			#print STDERR "\tAdded field 11, numFlds = " . scalar(@f) . "\n";
			
			
			
			my $hashStartRow = 7;
			my $hashEndRow = 8;
			my $rowKey = $pos;
			my $prior = $chrRowsH->{$rowKey};
			
			my $firstSeen = !$keysSeenThisFile->{$rowKey};
			$keysSeenThisFile->{$rowKey}++;
			

			# Dump if necessary
			if ($dumpRows)
			{
				if ($prior)
				{
					for (my $i=0; $i<=($hashEndRow-$hashStartRow); $i++)
					{
						#print STDERR "Dumping $rowKey .. Incrementing col " . ($i+$hashStartRow) . ": " . @{$prior}[$i] . "\n";
						$f[$i+$hashStartRow] = @{$prior}[$i];
					}
				}
				print W join("\t", @f);
				print W "\n";
			}	
			else
			{
				# Now we have to see if it has to be merged with another one, then store
				# it. (ONLY DO $firstSeen, i.e. the first one for each file)
				if ($firstSeen)
				{
					if ($prior)
					{
						for (my $i=0; $i<=($hashEndRow-$hashStartRow); $i++)
						{
							#print STDERR "Incrementing col " . ($i+$hashStartRow) . ": " . $prior . "\n";
							$f[$i+$hashStartRow] += @{$prior}[$i];
						}
					}
					#print STDERR "Storing entry for $rowKey: " . join(",",@f[$hashStartRow..$hashEndRow]) . "\n";
					my @smallF = @f[$hashStartRow..$hashEndRow];
					$chrRowsH->{$rowKey} = \@smallF;
				}
			}

		}
		
		close(R);
		$fileNum++;
	}
	
	close(W) if ($dumpRows);
	return $out;
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
	LINE: while (my $line = <F>)
	{
		$cpgNum++;
		if (($cpgNum % 1000000)==0)
		{
			print STDERR "On Cpg $cpgNum (chr=$lastChr, pos=$lastPos)\n";
		}
		
		## TESTING ## last LINE if ($cpgNum>100000); 
		
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



