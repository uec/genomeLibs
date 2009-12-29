#!/usr/bin/perl

use strict;

my $USAGE = "wrap_maq_mapmerge.pl out.map a.map b.map c.map ...";

# Argumetns
die "$USAGE\n" unless (@ARGV >= 2);
my $outfile = shift(@ARGV);
my @infiles = @ARGV;

# Check that input exists
my @goodInfiles = ();
foreach my $infile (@infiles)
{
	if (-s $infile)
	{
		push(@goodInfiles, $infile);
	}
	else
	{
    	print STDERR "wrap_maq_mapmerge: input file ${infile} is 0 length\n" ;
	}
}

# Run
my $cmd = join(" ", "maq","mapmerge",$outfile, @goodInfiles);
print STDERR "${cmd}\n";
print STDERR `${cmd} 2>&1`;

# Check that the output exists
die "wrap_maq_mapmerge: output file ${outfile} is 0 length\n" unless (-s $outfile);

