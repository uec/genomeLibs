#!/usr/bin/perl

use strict;
use File::Basename;

my @files = @ARGV;


foreach my $f (@files)
{
    my ($name, $path, $suf) = fileparse($f, qr/\.[^.]*/);
    my $outfn = join("", $path, $name, ".NOBADCHROMS", $suf);
    print "($name)\t($path)\t($suf)\n";

    die "Couldn't read $f\n" unless (open(IN,$f));
    die "Couldn't write to $outfn\n" unless (open(OUT, ">$outfn"));
    my $count = 1;
    while (my $line = <IN>)
    {
		chomp $line;	
	
		if (($line =~ /chr\d+(.)/) ||
			($line =~ /chrM(.)/) ||
			($line =~ /chrX(.)/) || 
			($line =~ /chrY(.)/))
		{
			my $nextChar = $1;
		    print OUT $line."\n" unless ($nextChar eq "_");
		} 
    }
    
    # Close input file
    close(IN);
    close(OUT);
}
