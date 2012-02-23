#!/usr/bin/env perl

use strict;
use Getopt::Long;

my $USAGE = "cat in.gtf | gtfFilterOutOffChromEntries hg18 > out.gtf";
my $GENOMES = "/home/uec-00/shared/production/genomes";
my $REF_FN_MAP = { "hg18" => "${GENOMES}/hg18_unmasked/hg18_unmasked.plusContam.dict",
	"hg19" => "${GENOMES}/hg19_rCRSchrm/hg19_rCRSchrm.dict" };

# my $distUpstream = 1000;
# my $minMapq = 20;
# GetOptions ('distUpstream=i', \$distUpstream, 'minq=i'=>\$minMapq) || die "$USAGE\n";

# Input params
die "$USAGE\n" unless (@ARGV==1);
my ($genome) = @ARGV;

# Get reference.
my $refFn = getRefFn($genome);

# Read lengths
my $lensByChrom = getChromLengths($refFn);

# Filter
LINE:while (my $line = <STDIN>)
{
    chomp $line;
    my @flds = split(/\t/,$line);

    my $lineGood = 1;
    if (scalar(@flds)>=8)
    {
	my ($chr, $source, $type, $s, $e, $score, $strand, $dummya, $dummyb) = @flds;
	$lineGood = 0 if ($s < 1);

	# Check end
	my $chromLen = $lensByChrom->{$chr};
	if (!$chromLen)
	{
	    print STDERR "Can't find length of chrom ${chr}!  Quitting\n";
	    die;
	}
	$lineGood = 0 if ($e >= ($chromLen-1));
    }

    print join("\t",@flds)."\n" if ($lineGood);
}

# - - - - Functions

sub getRefFn
{
	my ($genome) = @_;

# 	if ($ingtf =~ /hg18/i)
# 	{
# 		$genome = "hg18";
# 	}
# 	elsif ($ingtf =~ /hg19/i)
# 	{
# 		$genome = "hg19";
# 	}

	my $out = $REF_FN_MAP->{$genome};

	unless ($out)
 	{
 		print STDERR "Can't find reference dict file for genome \"$genome\"\n";
 		die;
 	}

	return $out;
	
}


sub getChromLengths
{
    my ($dictFn) = @_;

    my $out = {};
    die "Can't read dict file ${dictFn}\n" unless (open(F,$dictFn));
    while (my $line = <F>)
    {
	chomp $line;
	if ($line =~ /SN\:(\w+)\s+LN\:(\w+)\s/)
	{
	    my ($chr, $len) = ($1, $2);
#	    print STDERR "Chr=$chr\tlen=$len\n";
	    $out->{$chr} = $len;
	}
    }
    close(F);

    return $out;
}


