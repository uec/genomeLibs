#!/usr/bin/env perl

use strict;
use File::Basename;
use Getopt::Long;

my $USAGE = "maqToSeqTagsMysqlImp.pl --outPrefix mnase < input.txt";

my $delim = "\t";
my $outPrefix = "features";
my $desc = "";
my $h0unique = 0;
GetOptions ('outPrefix=s' => \$outPrefix, 'h0unique' => \$h0unique) || die "$USAGE\n";

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
	print STDERR "On line $lineNum\n" if (($lineNum%100000)==0);
	
	chomp $line;
	my @f = split(/\t/, $line);
	my $nF = scalar(@f);
	
	# Check for valid line
	if ($nF < 14)
	{
		print STDERR "Unrecognized SAM format! $nF fields\n$line\n";
		next LINE;
	} 
	
	if ($h0unique)
	{
		$line !~ /h0\:i\:(\d+)/i;
		my $h0 = $1;
		$line !~ /h1\:i\:(\d+)/i;
		my $h1 = $1;
		my $total = $h0+$h1;
		
		if ($total == 0)
		{
#			print "No hits: $line\n";
			next LINE;
		}
		elsif ($h0>0)
		{
			if ($h0>1)
			{
#				print "Non-unique seq: $line\n";
				next LINE;
			}
		}
		elsif ($h1>1)
		{
#			print "Non-unique seq: $line\n";
			next LINE;
		}
	}
	
	my $chr = $f[2];
	my $pos = $f[3];
	my $strand = ($f[1] & 16) ? "-" : "+";
	
	# THe pos is left-most, so we have to adjust for strand
	if ($strand eq '-')
	{
		my $seq = $f[9];
		my $len = length($seq);
		$pos -= $len;
	}
	
	
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

