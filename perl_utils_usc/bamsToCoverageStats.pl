#!/usr/bin/perl

use strict;
use File::Basename;
use Getopt::Long;

my $USAGE = "bamsToCoverageStats.pl chr1.bam chr2.bam ... > stats.csv";


#GetOptions ('chromCol=i' => \$chromCol, 'omitChromFldInOutput' => \$omitChrFld, 'outPrefix=s' => \$outPrefix, 'delimx=s' => \$delim) || die "$USAGE\n";


die "$USAGE\n" unless (@ARGV>0);


my $counts = {};
my $totalUnique = 0;
my $totalBases = 0;
foreach my $fn (@ARGV)
{
	my $cmd = "samtools view -ub -q 30 ./tumor_61AY8AAXX_s2345678.NODUPS.sorted.calmd.q20.IGF2_80kb.bam | samtools pileup -s - ";
	die "Couldn't pipe $cmd\n" unless (open(PILE, "${cmd} |"));
	
	while (my $line = <PILE>)
	{
		chomp $line;
		my @f = split(/\t/,$line);
		
		my $reads = $f[3];
		$counts->{$reads}++;
		
		$totalUnique++;
		$totalBases += $reads;
	} 
	close PILE;

	print STDERR "File ${fn}: ${totalUnique} unique, ${totalBases} total\n";
}

print STDOUT "${totalUnique} unique, ${totalBases} total\n";
foreach my $count (sort {$a<=>$b} (keys(%$counts)))
{
	print STDOUT "${count}," . $counts->{$count} . "\n";
}
 