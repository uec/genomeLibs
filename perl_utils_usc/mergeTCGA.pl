#!/usr/bin/env perl

use strict;
use Getopt::Long;
use File::Basename qw/basename/;

my $USAGE = "mergeTCGA.pl sampleList.csv aliquotFile.csv";

my $sampleListCol = 1;
my $aliquotFileCol = 1;
my $title = "subset";
my $dtypeonly = 0;
my $aliquotDelim = "\t";
GetOptions ('aliquotDelim=s',\$aliquotDelim,'sampleListCol=i', \$sampleListCol, 'aliquotFileCol=i',\$aliquotFileCol, 'title=s', \$title, 'dtypeonly!', \$dtypeonly) || die "$USAGE\n";

# Input params
die "$USAGE\n" unless (@ARGV==2);
my ($sampFn, $alFn) = @ARGV;

# get samples
my $samples = readSamples($sampFn,$sampleListCol);

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
    my $rawsamp = $f[$aliquotFileCol-1];
    my $outsamp = idToSample($rawsamp);

    my $ncols = scalar(@f);

    if (!$seenheader && (($line =~ /case/i) || ($line =~ /sample/i)))
    {
#	print STDERR "SAW HEADER: $line\n";
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
    }
}
close(INF);
close(OUTF);



# Put the header line in one called "header"
sub readSamples
{
    my ($fn, $col) = @_;
    die "Can't read $fn\n" unless (open(INF,$fn));
    my $samps = {};
    while (my $line = <INF>)
    {
	chomp $line;
	my @lines_fixed = split(chr(0xd), $line); # Remove windows linefeeds 
	$line = $lines_fixed[0];
	my @f = split(/,/,$line);
	@f = map {$a=$_; chomp $a; $a;} @f; # Remove windows linefeeds
	my $samp = $f[$col-1];
	if (!$samps->{header} && ($samp !~ /^\s+$/) && ($samp !~ /tcga/i))
	{
	    # Hopefully this is the header
	    print STDERR "Found 2nd set header\n";
	    $samp = "header";
	}
	else
	{
	    $samp = idToSample($samp);
#	    print "Found samp $samp\n";
	}
	$samps->{$samp} = join(",",@f);
    }
    close(INF);

    return $samps;
}

sub idToSample
{
    my ($id) = @_;
    my $sample = $id;

    if ($sample =~ /(TCGA-\w\w+-\w\w\w\w+)/)
    {
	$sample = $1;
    }

    return $sample;
}


