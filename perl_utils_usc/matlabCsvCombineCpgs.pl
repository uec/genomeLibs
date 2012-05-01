#!/usr/bin/env perl

use File::Basename;
use File::Temp qw/ tempfile tempdir /;
use strict;
use Getopt::Long;

my $USAGE = "matlabCsvCombineCpgs.pl file1.csv file2.csv file3.csv ...";
GetOptions () || die "$USAGE\n";
# 'compareCutoff=f',\$compareCutoff,'testCutoff=f',\$testCutoff, 
#print STDERR "doBare=${doBare}\tintervalFile=${intervalFile}\n";

my $IND_CHR = 0;
my $IND_POS = 1;
my $IND_STRAND = 2;
my $IND_METH = 3;
my $IND_COUNT = 4;

# Input params
die "$USAGE\n" unless (@ARGV>=1);
my (@fns) = @ARGV;


foreach my $fn (@fns)
{
    my $base = basename($fn, qw/.csv/);
    my $fnout = "${base}.combinedCG.csv";

    print STDERR "$fn -> $fnout\n";

    die "Can't read $fn\n" unless (open(FIN,$fn));
    die "Can't write to $fnout\n" unless (open(FOUT,">$fnout"));

    convertFile(\*FIN,\*FOUT);

    close(FOUT);
    close(FIN);
}


sub convertFile
{
    my ($fin, $fout) = @_;

    my $fPrev = [];
    my $pairsFound = 0;
    while (my $line=<$fin>)
    {
	chomp $line;
	my $fNext = [split(/,/,$line)];
	changeToFw($fNext);
	print STDERR "Following line doesn't have 5 fields: $line\n" unless (scalar(@$fNext)==5);

	# Always output the "previous" one.  Unless it's a pair, in which case we replace the sum.
	if (cpgPair($fPrev,$fNext))
	{
	    # Found a pair
	    $pairsFound++;

	    my @fCombined = combinePair($fPrev, $fNext);
	    #print STDERR "Found pair, a=(".join(",",@$fPrev).")\tb=(".join(",",@$fNext).")\t-->\tcombined=(".join(",",@fCombined).")\n";

	    $fNext = \@fCombined;
	}
	else
	{
	    # Not a pair.  Output
	    print $fout join(",",@$fPrev)."\n";
	}

	# And replace prev
	$fPrev = $fNext;
    }

    # We still have to output the final one
    print $fout join(",",@$fPrev)."\n";
    
    print STDERR "Found ${pairsFound} proper pairs\n";
}

sub combinePair
{
    my ($fPrv, $fNxt) = @_;
    
    my @out = @$fPrv; # Make everything match the forward strand one
    #@out[$IND_STRAND] = ".";
    
    my $totalM = numMeth($fPrv) + numMeth($fNxt);
    my $total =  numReads($fPrv) + numReads($fNxt);

    my $meth = sprintf("%.1f",($totalM/$total)*100);
    
    @out[$IND_METH] = $meth;
    @out[$IND_COUNT] = $total;

    return @out;
}

sub changeToFw($f)
{
    my ($f) = @_;

    if (@{$f}[$IND_STRAND]==-1)
    {
	@$f[$IND_STRAND] = 1;
	@$f[$IND_POS]--;
    }
}

sub numMeth
{
    my ($f) = @_;

    my $nm = sprintf("%.1f",@{$f}[$IND_METH] * @{$f}[$IND_COUNT] / 100.0);
    return $nm;
}

sub numReads
{
    my ($f) = @_;

    my $n = @{$f}[$IND_COUNT];
    return $n;
}

# Assumes the neg strand has already been converted to + 
sub cpgPair
{
    my ($fPrv, $fNxt) = @_;

    my $paired=  ( (@{$fPrv}[$IND_CHR] == @{$fNxt}[$IND_CHR]) &&
		   (@{$fPrv}[$IND_POS] == @{$fNxt}[$IND_POS]) );

    return $paired;
}


sub cpgPairStranded
{
    my ($fPrv, $fNxt) = @_;

    my $paired=  ( (@{$fPrv}[$IND_CHR] == @{$fNxt}[$IND_CHR]) &&
		   (@{$fPrv}[$IND_POS] == (@{$fNxt}[$IND_POS] - 1)) &&
		   (@{$fPrv}[$IND_STRAND] == 1) &&
		   (@{$fNxt}[$IND_STRAND] == -1) );

    return $paired;
}
