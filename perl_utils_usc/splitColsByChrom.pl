#!/usr/bin/perl

use strict;
use File::Basename;
use Getopt::Long;

my $USAGE = "splitColsByChrom.pl --omitChromFldInOutput --chromCol 1 --delim , --outPrefix fileout file1.txt file2.csv\nchromCol numbering starts at 1";

my $delim = "\t";
my $outPrefix = 0;
my $chromCol = 1;
my $omitChrFld = 0;
GetOptions ('chromCol=i' => \$chromCol, 'omitChromFldInOutput' => \$omitChrFld, 'outPrefix=s' => \$outPrefix, 'delimx=s' => \$delim) || die "$USAGE\n";


die "$USAGE\n" unless (@ARGV>0);
foreach my $fn (@ARGV)
{
    my ($name, $path, $suf) = fileparse($fn, qr/\.[^.]*/);

	$outPrefix = $name unless ($outPrefix);
	print STDERR "chromCol=$chromCol\tprefix=$outPrefix\tdelim=\"$delim\"\tfile=$name\n";
	
	die "Can't read $fn\n" unless (open(F,$fn));
	my $outfhByChrom = {};
	
	
	while (my $line = <F>)
	{
		chomp $line;
		my @flds = split($delim, $line);
		
		my $chr = $flds[$chromCol-1];
		$chr = "chr${chr}" unless ($chr =~ /^chr/i);
		$flds[$chromCol-1] = $chr;
		
		splice (@flds, $chromCol-1, 1) if ($omitChrFld);
		
		my $fh = $outfhByChrom->{$chr};
		if ($fh)
		{
			#print STDERR "Found fh for chr $chr: $fh\n";	
		}
		else
		{
			#print STDERR "Didn't find fh for chr $chr\n";	
			my $outfn = "${outPrefix}${chr}.txt";
			die "Can't write to $outfn\n" unless (open($fh, ">$outfn"));
			$outfhByChrom->{$chr} = $fh;
		}
		
		print $fh join($delim, @flds)."\n";
	}
	close(F);

	
	# Close output files	
	foreach my $outfh (values(%$outfhByChrom))
	{
		close($outfh);
	}
	

}