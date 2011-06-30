#!/usr/bin/perl

use strict;
use File::Basename;

my $USAGE = "gffPullRelativeRegions.pl -100 40 file1.gtf file2.gtf ...";

die "Must have at least 3 command line params:\n$USAGE\n" unless (scalar(@ARGV)>=3);
my ($startRel, $endRel, @files) = @ARGV;

die "end must be greater than start:\n$USAGE\n" unless ($endRel>=$startRel);

foreach my $f (@files)
{

    my ($name, $path, $suf) = fileparse($f, qr/\.[^.]*/);
    my $outfn = sprintf("%s%s.RELPOS%sTO%s%s",$path,$name,addsign($startRel),addsign($endRel),$suf);
    print STDERR "$f --> $outfn\n";

    die "Couldn't read $f\n" unless (open(IN,$f));
    die "Couldn't write to $outfn\n" unless (open(OUT, ">$outfn"));
    my $count = 1;
    while (my $line = <IN>)
    {
	chomp $line;
	if ($line =~ /^\s*\#/ || $line =~ /^\s*$/ || $line =~ /track/)
	{
	    print OUT $line."\n";
	}
	elsif ($line =~ /random/i)
	{
	    # Remove chromosome random
	}
	else
	{
	    my @flds = split(/\t/,$line);

	    # Adjust start and end fields
	    my $rev = ($flds[6] eq "-");
	    my $midpoint = int(($flds[3]+$flds[4])/2);
	    my $news = $midpoint + (($rev) ? (-$endRel) : $startRel);
	    my $newe = $midpoint + (($rev) ? (-$startRel) : $endRel);
	    $flds[3] = $news;
	    $flds[4] = $newe;
	     
	    # Fix ID field.
	    my $id = @flds[8];
	    
	    my $new_id = $count++;
	    if ($id =~ /\w+\s+\"([^\"]+)\"\;/)
	    {
		$new_id .= "-".$1;
	    }
	    @flds[8] = "gene_id \"$new_id\"; transcript_id \"$new_id\";";


	    # Fix chr field if necessary
	    @flds[0] =~ s/chr23/chrX/g;
	    @flds[0] =~ s/chr24/chrY/g;
	    

	    # We also see if the source field needs correcting
	    if (@flds[1] =~ /^[\d\,]+$/)
	    {
		@flds[1] = "gff_fix_id";
	    }

	    # Also sometimes the 3rd field gets screwed up
	    @flds[2] = "exon";

	    # and print
	    print OUT join("\t",@flds)."\n";
	}
    }
    
    # Close input file
    close(IN);
    close(OUT);
}

sub addsign
{
    my ($num) = @_;

    my $out = $num;
    if ($num>=0)
    {
	$out = "+${num}";
    }
    return $out;
}

