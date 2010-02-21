#!/usr/bin/env perl

use strict;
use Getopt::Long;

my $USAGE = "wigScanPeaks.pl -minScore 0.05 -minCount 1 file1.wig file2.wig ...";

my $minScore = -1;
my $minCount = 1;
my $span = 100;
GetOptions ('minScore=f' => \$minScore, 'minCount=i' => \$minCount, 'span=i' => \$span, ) || die "$USAGE\n";
die "$USAGE\n" unless (($minScore>=0) && (@ARGV>0));

foreach my $fn (@ARGV)
{
    my $base = $fn;
    $base =~ s/\.wig.*//g;
    my $scoresec = ".minScore${minScore}";
    my $countsec = ".minCount${minCount}";
    my $spansec = ".span${span}";
    my $outfn = "${base}${scoresec}${countsec}${spansec}.gtf";
	
    print STDERR "$fn -> $outfn\n";
    die "Can't read $fn\n" unless (open(IN,$fn));
    die "Can't write to $outfn\n" unless (open(OUT,">$outfn"));
    
    my $halfspan = int($span/2);

    my $on_line = 0;
    my $curchr = 0;
    my $winds = -1;
    my $winde = -1;
    my $windcount = 0;
    LINE: while (my $line = <IN>)
    {
		chomp $line;
		
		if ($line =~ /chrom=(chr.*)/)
		{
			# Check if we have a previous window
			checkAndProcessWind(\*OUT, $winds, $winde, $windcount, $minCount, $curchr);
			
			
			$curchr = $1;
			print STDERR "New chrom: $curchr\n";
			$winds = -1;
			$winde = -1;
			$windcount = 0;
		}
		else
		{
			my @flds = split(/\t/,$line);
			if (scalar(@flds)!=2)
			{
				print STDERR "Illegal line: $line\n";
			}
			else
			{
				my ($mid, $score) = @flds;
				my $s = $mid-$halfspan;
				my $e = $mid+$halfspan;
				
				my $pass = ($score >= $minScore);
				
				# If we're not in a window, process our current window
				if (!$pass)
				{
					checkAndProcessWind(\*OUT, $winds, $winde, $windcount, $minCount, $curchr);
					$winds = -1;
					$winde = -1;
					$windcount = 0;
				}
				else
				{
					# We're in a window.  Update all vars
					$winds = $s if ($winds==-1);
					$winde = $e;
					$windcount++;
				}
			}
		}		

	
		$on_line++;
    }

	# See if we have a final window
	checkAndProcessWind(\*OUT, $winds, $winde, $windcount, $minCount, $curchr);
	

    close(OUT);
    close(IN);
}

sub checkAndProcessWind
{
	#print STDERR join(", ", @_)."\n";
	my ($outfh, $winds, $winde, $windcount, $minCount, $curchr) = @_;
	
	# Check that it's valid
	return if ($winds == -1);
	return if ($windcount < $minCount);
	
	$winds = 1 if ($winds==0); # 1-based
	
	my $id = "${curchr}:${winds}-${winde}";
	
	# Write wind
	print $outfh join("\t", 
		$curchr, 
		"wigscanpeaks",
		"exon",
		$winds,
		$winde,
		$windcount,
		".",
		".",
		"gene_id \"$id\" ;transcript_id \"$id\"") . "\n";
}
