#!/usr/bin/env perl

use strict;
use Getopt::Long;
use File::Basename qw/basename/;

my $USAGE = "mergeTCGA.pl sampleList.csv aliquotFile.csv -aliquotDelim , -sampleListCol 1 -aliquotFileCol 1 -title WGBS -dtypeonly -useDnaType";

my $sampleListCol = 1;
my $aliquotFileCol = 1;
my $title = "subset";
my $dtypeonly = 0;
my $aliquotDelim = "\t";
my $useDnaType = 0;
GetOptions ('aliquotDelim=s',\$aliquotDelim,'sampleListCol=i', \$sampleListCol, 'aliquotFileCol=i',\$aliquotFileCol, 
		'title=s', \$title, 'dtypeonly!', \$dtypeonly, 'useDnaType!', \$useDnaType) || die "$USAGE\n";

# Input params
die "$USAGE\n" unless (@ARGV==2);
my ($sampFn, $alFn) = @ARGV;

# get samples
my $samples = readSamples($sampFn,$sampleListCol, $aliquotDelim, $useDnaType);

# Now output subset
my $outfn = basename($alFn);
$outfn =~ s/\.csv//g;
$outfn .= ".${title}.csv";

die "Can't write to $outfn\n" unless (open(OUTF,">$outfn"));
die "Can't read from $alFn\n" unless (open(INF,$alFn));
my $seenheader = 0;
my $maxcols = -1;
while (my $line = <INF>)
{
    chomp $line;
    my @lines_fixed = split(chr(0xd), $line); # Remove windows linefeeds 
    $line = $lines_fixed[0];

    my @f = split(/${aliquotDelim}/,$line);
    @f = split(/,/,$line) if (scalar(@f)==1); # Special catch
    
    my $rawsamp = $f[$aliquotFileCol-1];
    my $outsamp = idToSample($rawsamp, $useDnaType);
    # Chek if we want the non DNA type one
    if ($useDnaType && !$samples->{$outsamp})
    {
    	$outsamp = idToSample($rawsamp, 0);
    }

    my $ncols = scalar(@f);
#    print STDERR sprintf("Found %d cols\n",$ncols);

    if (!$seenheader && (($line =~ /case/i) || ($line =~ /sample/i)))
    {
	print STDERR "SAW HEADER: $line\n";
	# Header
	print OUTF join(",",@f,$samples->{"header"})."\n";
	$maxcols = $ncols;
	$seenheader = 1;
    }
    elsif ($dtypeonly && ($rawsamp !~ /D$/i))
    {
    }
    else
    {
#	print STDERR "Checking samp ${rawsamp}/${outsamp}\n";
		if ($samples->{$outsamp})
		{
		    print OUTF join(",",map { $a=$_; $a=~s/\"/__/g; ($a=~/\,/) ? "\"$a\"" : "$a"} @f[0..$maxcols-1]).",";
		    print OUTF join(",", map {$_} 1..($maxcols-$ncols))."," if ($maxcols>$ncols);
		    print OUTF join(",",$samples->{$outsamp});
		    print OUTF  "\n";
		}
		else
		{	
			print STDERR "Couldn't find sample $outsamp\n";
		}
    }
}
close(INF);
close(OUTF);



# Put the header line in one called "header"
sub readSamples
{
    my ($fn, $col, $delim, $useDnaType) = @_;
    
    $delim = "," unless ($delim);
    
    die "Can't read $fn\n" unless (open(INF,$fn));
    my $samps = {};
    while (my $line = <INF>)
    {
	chomp $line;
	my @lines_fixed = split(chr(0xd), $line); # Remove windows linefeeds 
	$line = $lines_fixed[0];
	
	my @f = split(/${delim}/,$line);
    @f = split(/,/,$line) if (scalar(@f)==1); # Special catch

	@f = map {$a=$_; chomp $a; $a;} @f; # Remove windows linefeeds
	my $samp = $f[$col-1];
	if (!$samps->{header} && ($samp !~ /^\s+$/) && ($samp !~ /tcga\-/i))
	{
	    # Hopefully this is the header
	    print STDERR "Found 2nd set header\n";
	    $samp = "header";
	}
	else
	{
	    $samp = idToSample($samp, $useDnaType);
#	    print "Found samp $samp\n";
	}
	$samps->{$samp} = join(",",@f);
    }
    close(INF);

    return $samps;
}

sub idToSample
{
    my ($id, $useDnaType) = @_;
    my $sample = $id;

	# Special cases
	$sample =~ s/^WU\-GBM/TCGA/g;
	$sample =~ s/^GBM/TCGA/g;

    my $outsample = 0;
	if ($useDnaType && ($sample =~ /(TCGA-\w\w+-\w\w\w\w+-\w\w\w+)/))
	{
		$outsample = $1;
	}
	elsif ($sample =~ /(TCGA-\w\w+-\w\w\w\w+)/)
	    {
			$outsample = $1;
	    }
	
	print STDERR sprintf("%s --> %s\n", $id, $outsample);
    return $outsample;
}


