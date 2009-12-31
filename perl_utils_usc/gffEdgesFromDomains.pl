#!/usr/bin/env perl

use strict;
use File::Basename;
use Getopt::Long;

my $USAGE = "gffEdgesFromDomains.pl --useOnlyFivePrimeEnd in1.gtf in2.gtf ... ";

my $useOnlyFivePrimeEnd = 0;
GetOptions ('useOnlyFivePrimeEnd!' => \$useOnlyFivePrimeEnd) || die "$USAGE\n";


die "$USAGE\n" unless (@ARGV>0);
foreach my $fn (@ARGV)
{
	die "Can't read $fn\n" unless (open(INF,$fn));

    my ($name, $path, $suf) = fileparse($fn, qr/\.[^.]*/);
    my $filelabel = ($useOnlyFivePrimeEnd) ? "StartEdges" : "StartEndEdges";
	my $outfn = "${name}.${filelabel}.gtf";
	die "Can't write to $outfn\n" unless (open(OUTF,">$outfn"));

	while (my $line=<INF>)
	{
	    chomp $line;
	    my @flds = split(/\t/,$line);
	
	    if (@flds<8)
	    {
	        print STDERR "Illegal: $line\n";
	        print OUTF "$line\n";
	    }
	    else
	    {
	        my $rangeS = $flds[3];
	        my $rangeE = $flds[4];
	        
	        my $instrand = $flds[6];

	        for (my $i = 0; $i <= 1; $i++)
	        {
				# i = 0 => 5', i=1 => 3'
				my $label = ($i==0) ? "5prime" : "3prime";
				
				my @newflds = @flds; # Make a copy
				$newflds[6] = ($i==0) ? "+" : "-";
				$newflds[3] = ($i==0) ? $rangeS : $rangeE;	
				$newflds[4] = ($i==0) ? $rangeS : $rangeE;	
				$newflds[8] =~ s/\"([^\"]+)\"/\"$1-$label\"/g if ($flds[8]);
				
				# Do we want to actually include this one?
				my $correctStrand = ($instrand eq (($i==0) ? "+" : "-"));
				#print STDERR "correct=${correctStrand}\n";
				if (!$useOnlyFivePrimeEnd || ($instrand eq ".") || $correctStrand)
				{
					print OUTF join("\t", @newflds)."\n";	
				}
				
	        }  
	    }      
	}

	close (OUTF);
	close (INF);
}
