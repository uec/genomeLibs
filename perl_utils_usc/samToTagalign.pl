#!/usr/bin/env perl

use strict;
use File::Basename;
use Getopt::Long;

my $USAGE = "cat input.sam | samToBed.pl > out.tagAlign";
# GetOptions ('outPrefix=s' => \$outPrefix, 'h0unique' => \$h0unique, 'noExtras' => \$::noExtras) || die "$USAGE\n";


die "$USAGE\n" unless (@ARGV==0);


my $lineNum = 0;
LINE: while (my $line = <STDIN>)
{
	$lineNum++;
	print STDERR "On line $lineNum\n" if (($lineNum%100000)==0);

	next if ($line =~ /^\@/); # Header
	
	chomp $line;
	my @f = split(/\t/, $line);
	my $nF = scalar(@f);
	
	# Check for valid line
	if ($nF < 11)
	{
		print STDERR "Unrecognized SAM format! $nF fields\n$line\n";
		next LINE;
	} 
	
	my $chr = $f[2];
	my $pos = $f[3];
	my $strand = ($f[1] & 16) ? "-" : "+";

	my $seq = $f[9];
	my $len = length($seq);
	my $mapq = $f[4];


	print STDOUT join("\t",
			  $chr,
			  $pos,
			  $pos+$len-1,
			  $seq,
			  int(1000*$mapq/40),
			  $strand);
	print STDOUT "\n";
}


