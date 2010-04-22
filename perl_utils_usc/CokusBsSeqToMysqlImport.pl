#!/usr/bin/env perl

use strict;
use File::Basename;
use Getopt::Long;

my $USAGE = "CokusBsSeqToMysqlImport.pl --outPrefix mnase < input.txt";

my $delim = "\t";
my $outPrefix = "features";
my $desc = "";
GetOptions ('outPrefix=s' => \$outPrefix) || die "$USAGE\n";

die "$USAGE\n" unless (@ARGV==0);

my $CPG_COORDS_PREFIX = "/Users/benb/genomic-data-misc/CpGs/mm9/";
my $CPG_COORDS_SUFFIX = ".fa.cpgCoords.csv";


# Setup output files
$::openFhs = {};

my $fieldmap = {};
my $lineNum = 0;
my $seen = {};
my $numDups = 0;
my $lastKey = "";
my $lastChr = 0;
my $curTotal = 0;
my $curMeth = 0;
my $cgMap = {};
my $methylNonCgs = 0;
LINE: while (my $line = <STDIN>)
{
	$lineNum++;
	print STDERR "On line $lineNum\n" if (($lineNum%100000)==0);
	
	chomp $line;
	my @f = split(/\t/, $line);
	
	# Check for valid line
	if (@f != 2)
	{
		print STDERR "Unrecognized Cokus format! $line\n";
		next LINE;
	} 
	
	my $m = lc($f[1]);
	my ($chr, $pos, $strand);
	if ($f[0] =~ /m?(\d+)([\-\+])(\d+)/)
	{
		$chr = "chr" . int($1);
		($strand,$pos) = ($2,int($3));
	}
	elsif ($f[0] =~ /m?([XYMxym])[XYMxym]([\-\+])(\d+)/)
	{
		$chr = "chr${1}";
		($strand,$pos) = ($2,int($3));
	}
	elsif ($f[0] =~ /unk/ || $f[0] =~ /rnd/ )
	{
		next LINE;
	}
	else
	{
		print STDERR "Unrecognized Cokus format! $line\n";
		next LINE;
	}
	

	# Check if it's the last one
	my $key = join("__",$chr, $pos,$strand);
	print STDERR "Key: ${key}\n" unless ($lastChr eq $chr);
	
	# If it's a new chromsome, read in CpGs
	if ($chr ne $lastChr)
	{
		$cgMap = readCgMap($chr);
	}
	
	
	# Is it a CpG?
	
	if (!$cgMap->{$pos})
	{
		# Doesn't map.  Count methylated nonmaps
		$methylNonCgs++ if ($m eq 'm');
	}
	else
	{
		
		if ($lastKey eq $key)
		{
	#		print STDERR "\tFound a dup: $key\n";
		}
		else
		{
			# End the old one
			my ($lastChr, $lastPos, $lastStrand) = split(/__/,$lastKey);
			writeRow($outPrefix, $lastChr, $lastPos, $lastStrand, $curTotal, $curMeth);
			
			# start a new one
	#		print STDERR "\tFound a new one: $key\n";
			$curTotal = 0;
			$curMeth = 0;
		}
	
		# Increment
		$curTotal++;
		$curMeth++ if ($m eq 'm');
	}
	


	$lastChr = $chr;
	$lastKey = $key;
}

print STDERR "Saw ${numDups} duplicates in ${lineNum} reads\n";
print STDERR "Saw ${methylNonCgs} non-CG methyls\n";

# close files
foreach my $fh (values(%$::openFhs))
{
	close($fh);
}

# - - - - Functions

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

sub readCgMap
{
	my ($chr) = @_;
	
	my $fn = $CPG_COORDS_PREFIX.$chr.$CPG_COORDS_SUFFIX;
	print STDERR "Reading CGs from file $fn\n";
	
	die "Can't read $fn\n" unless (open(F,$fn));
	
	my $map = {};
	while (my $line = <F>)
	{
		chomp $line;
		my ($fChr, $pos) = split(/,/,$line);
		$map->{$pos}++; 
	}
		
	close(F);
	
	return $map;
}

