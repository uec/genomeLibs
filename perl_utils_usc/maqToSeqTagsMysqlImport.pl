#!/usr/bin/env perl

use strict;
use File::Basename;
use Getopt::Long;

my $USAGE = "maqToSeqTagsMysqlImp.pl --outPrefix mnase < input.txt";

my $delim = "\t";
my $outPrefix = "features";
my $desc = "";
GetOptions ('outPrefix=s' => \$outPrefix) || die "$USAGE\n";

die "$USAGE\n" unless (@ARGV==0);


# Setup output files
$::openFhs = {};

my $fieldmap = {};
my $lineNum = 0;
my $seen = {};
my $numDups = 0;
LINE: while (my $line = <STDIN>)
{
	$lineNum++;
	print STDERR "On line $lineNum\n" if (($lineNum%1000)==0);
	
	chomp $line;
	my @f = split(/\t/, $line);
	
	# Check for valid line
	if (@f != 16)
	{
		print STDERR "Unrecognized Maq format!\n";
		next LINE;
	} 
	
	my $chr = $f[1];
	my $pos = $f[2];
	my $strand = $f[3];
	
	# Remove dups
	my $key = join("__",$chr, $pos,$strand);
	if ($seen->{$key})
	{
		$numDups++;
	}
	else
	{
		#output
		writeRow($outPrefix, $chr, $pos, $strand);
	}
	$seen->{$key}++;
}

print STDERR "Saw ${numDups} duplicates in ${lineNum} reads\n";

# close files
foreach my $fh (values(%$::openFhs))
{
	close($fh);
}

# - - - - Functions

sub writeRow
{
	my ($prefix, $chr, $pos, $strand) = @_;
	my $fn = "${prefix}_${chr}.txt";
	my $fh = $::openFhs->{$fn};
	if (!$fh)
	{
		die "Can't write to $fn\n" unless (open($fh,">$fn"));
		$::openFhs->{$fn} = $fh; 
	}
	
	print $fh join("\t",$pos, $strand, 1, 1, 0, 0, 0, 0, 0, 0, 0);
	print $fh "\n";
}

