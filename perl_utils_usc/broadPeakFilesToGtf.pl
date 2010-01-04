#!/usr/bin/env perl

use strict;
use File::Basename;
use Getopt::Long;

my $USAGE = "broadPeakFilesToGtf.pl in1.gtf in2.gtf ..";

#my $delim = "\t";
#my $outPrefix = 0;
#my $chromCol = 1;
#my $omitChrFld = 0;
#GetOptions ('chromCol=i' => \$chromCol, 'omitChromFldInOutput' => \$omitChrFld, 'outPrefix=s' => \$outPrefix, 'delimx=s' => \$delim) || die "$USAGE\n";


die "$USAGE\n" unless (@ARGV>0);
foreach my $fn (@ARGV)
{
    my ($name, $path, $suf) = fileparse($fn, qr/\.[^.]*/);

	my $outFn = $name . "${suf}.hg18.nodups.gtf";
	
	die "Can't read $fn\n" unless (open(F,$fn));
	die "Can't write to $outFn\n" unless (open(OUTF,">$outFn"));
	

	my $onPeak = 1;	
	while (my $line = <F>)
	{
		chomp $line;
		my @f = split(/\t/, $line);

		print OUTF join("\t",
			$f[0], "broadPeaks", "exon", $f[1], $f[2], $f[7], ".", ".", "gene_id \"peak${onPeak}\"; transcript_id \"peak${onPeak}\"") . "\n";

		$onPeak++;

	}
	close(OUTF);
	close(F);

	

}