#!/usr/bin/env perl

use strict;
use File::Basename;
use List::Util qw(min max);

# Args
my @files = @ARGV;


foreach my $f (@files)
{
    my ($name, $path, $suf) = fileparse($f, qr/\.[^.]*/);
    print "($name)\t($path)\t($suf)\n";

    # Open output file
    my $outfn = $path.$name;
    $outfn .= "_transformed" if ($suf =~ /bed/i);
    $outfn .= ".bed";
    print "Writing to $outfn\n";
    die "Can't write to $outfn\n" unless (open(OUT,">$outfn"));

    # Keep output strings in a hash
    my $linesByKey = {};

    die "Couldn't read $f\n" unless (open(F,$f));
    my $on_line = 1;
    while (my $line = <F>)
    {
	if ($line =~ /track/i)
	{
	    print OUT $line;
	    next;
	}

	my $key = lineToKey($line);
	$linesByKey->{$key} = $line;

	print STDERR "Finished line $on_line\n" if (($on_line % 100000)==0);
	$on_line++;
    }

    print STDERR "Finished reading $on_line lines. Sorting...\n";
    my $on_key = 1;
    foreach my $key (sort keys(%$linesByKey))
    {
	my $line = $linesByKey->{$key};
	print OUT $line;
	
	print STDERR "Finished writing $on_key\n" if (($on_key % 100000)==0);
	$on_line++;
    }

    close(OUT);
}

sub lineToKey
{
    my ($line) = @_;

    my @f = split(/\t/,$line);
    my $key = sprintf("%s-%10d-%10d",$f[0],$f[1],$f[2]);
    return $key;
}
