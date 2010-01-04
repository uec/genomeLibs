#!/usr/bin/perl

use strict;
use File::Basename;
use Getopt::Long;

my $USAGE = "gtfRemoveExactDups.pl in1.gtf in2.gtf ..";

#my $delim = "\t";
#my $outPrefix = 0;
#my $chromCol = 1;
#my $omitChrFld = 0;
#GetOptions ('chromCol=i' => \$chromCol, 'omitChromFldInOutput' => \$omitChrFld, 'outPrefix=s' => \$outPrefix, 'delimx=s' => \$delim) || die "$USAGE\n";


die "$USAGE\n" unless (@ARGV>0);
foreach my $fn (@ARGV)
{
    my ($name, $path, $suf) = fileparse($fn, qr/\.[^.]*/);

	my $outFn = $name . ".nodups.gtf";
	
	die "Can't read $fn\n" unless (open(F,$fn));
	die "Can't write to $outFn\n" unless (open(OUTF,">$outFn"));
	
	
	my $seen = {};
	while (my $line = <F>)
	{
		chomp $line;
		my @f = split(/\t/, $line);
		
		my $key = join("__",$f[0],$f[3],$f[4]);
		
		if ($seen->{$key})
		{
			#print STDERR "Already saw key $key\n";
		}
		else
		{
			print OUTF $line . "\n";
		}
		
		
		$seen->{$key}++;

	}
	close(OUTF);
	close(F);

	

}