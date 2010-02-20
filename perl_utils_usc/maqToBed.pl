#!/usr/bin/env perl

use strict;
use Getopt::Long;

my $USAGE = "maqToBed.pl -fragSize 400 (extends pos to halfway into frag) file1.map.q30.txt file2.map.q30.txt ...";

my $fragSize = 0;
my $noExtras = 0;
GetOptions ('fragSize=i' => \$fragSize, 'noextras' => \$noExtras) || die "$USAGE\n";
die "$USAGE\n" unless (@ARGV>0);

foreach my $fn (@ARGV)
{
    my $base = $fn;
    $base =~ s/\.map.*//g;
    my $fragsec = ($fragSize!=0) ? ".frag${fragSize}" : "";
    my $outfn = "${base}${fragsec}.bed";
	
	my $halfFrag = int($fragSize/2);

    print STDERR "$fn -> $outfn\n";
    die "Can't read $fn\n" unless (open(IN,$fn));
    die "Can't write to $outfn\n" unless (open(OUT,">$outfn"));

    print OUT "track name='${base}' visibility=4 itemRgb='on'\n" unless ($noExtras);
    my $on_line = 0;
    LINE: while (my $line = <IN>)
    {
		chomp $line;
		my @flds = split(/\t/,$line);
		my $name = $on_line; # $flds[0]
		my $chrom = $flds[1];
		next LINE unless ($chrom =~ /^chr/);
	
		my $s = $flds[2];
		my $strand = $flds[3];
		my $score = $flds[6];
		my $len = $flds[13];
	
		my $e = $s + $len - 1;
		if ($fragSize != 0)
		{
			# We place it at the middle of the estimated fragment.
			my $mid;
			if ($strand eq '+')
			{
				$mid = $s + $halfFrag;
			}
			else
			{
				$mid = $s + $len - $halfFrag;
				$mid = 1 if ($mid<1);
			}
			$s = $mid;
			$e = $mid;
		}
		
		
		my $color = ($strand eq '+') ? "0,0,255" : "255,0,0";
	
		my @out_flds = ( $chrom, $s, $e, $name, $score, $strand);
		push(@out_flds, $s, $e, $color) unless ($noExtras);
		print OUT join("\t",@out_flds)."\n";
	
		$on_line++;
    }

    close(OUT);
    close(IN);

    `bzip2 -f $outfn` unless ($noExtras);

}
