#!/usr/bin/perl

use strict;
use File::Basename;

my @files = @ARGV;

$::MIN_SCORE = 0;

foreach my $f (@files)
{

    my ($name, $path, $suf) = fileparse($f, qr/\.[^.]*/);
    print "($name)\t($path)\t($suf)\n";

    # Keep output strings in a hash
    my $src_h = {};

    die "Couldn't read $f\n" unless (open(F,$f));
    my $on_line = 1;
    while (my $line = <F>)
    {
	if ($line =~ /^\s*\#/ || $line =~ /^\s*$/ || $line =~ /track name/)
	{
	}
	else
	{
	    my @flds = split(/\t/,$line);
	    if (($flds[5] eq '.') || ($flds[5]>$::MIN_SCORE))
	    {
	    my $src = @flds[2];
	    $src =~ s/[\/]//g;
	    $src =~ s/\./-/g;
	    $src_h->{$src} .= $line;
	    }
	}
	
	print STDERR "On line $on_line\n" if (($on_line%10000)==0);
 	$on_line++;
    }
    
    # Close input file
    close(F);

    # And write output files
    SRC: foreach my $src (keys(%$src_h))
    {
	next SRC if (!$src || ($src =~ /^\s*\#/) || ($src =~ /^\s*$/));

	my $str = $src_h->{$src};
	$src =~ s/\/\#\://g;
	my $minSec = ($::MIN_SCORE>0) ? ".minScore".$::MIN_SCORE : "";
	my $outfn = $path.$name."_$src".$minSec.$suf;
	
	# Special cases
	$outfn =~ s/wgEncodeRegTfbsClustered_/ENCODEclst_/g;
	$outfn =~ s/minScore/min/g;
	$outfn =~ s/.gtf//g;
	
	
	
	print "Writing to $outfn\n";

	if (1)
	{
	    die "Can't write to $outfn\n" unless (open(OUT,">$outfn"));
	    print OUT $str;
	    close(OUT);
	}
    }
}
