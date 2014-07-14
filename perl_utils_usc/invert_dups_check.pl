#!/usr/bin/perl -w
##This script count the inverted dups by count the first N bases' match rate in two ends' sequences.

## take care, 1st end, T should be also considered as C, while in 2nd end, A should be also considered as G

## author: Yaping Liu  lyping1986@gmail.com 
## time: 2013-1-8

#Usage:  perl invert_dups_check.pl [option] output_log.txt sample.1st_end.fastq sample.2nd_end.fastq 


use strict;
use Getopt::Long;

my $output_log = $ARGV[0];
my $first_end = $ARGV[1];
my $second_end = $ARGV[2];

my $mismatch=3;
my $first_N=40;

sub usage {
	
    print "\nUsage:\n";
    print "perl invert_dups_check.pl [Options]  output_log.txt sample.1st_end.fastq sample.2nd_end.fastq\n\n";

    print " This script count the inverted dups by count the first N bases' match rate between two ends' sequences\n";
    print "It provides another way to count inverted dups other than count the same starting coordinate for mapped & proper paired reads, which could count inverted dups in all of the reads\n\n";
	print "[Options]:\n\n";
	
	print "  --mismatch INT : maximum mismathces allowed for account the inverted dups. (Default: 3)\n\n";
	print "  --first_N INT : First N bases to count for the matches. (Default: 40)\n\n";
}

GetOptions(
	"mismatch=i" => \$mismatch,
	"first_N=i" => \$first_N,
);

usage() if ( scalar(@ARGV) <= 3 );


my $linecount_within = 1;
my $linecount_global = 1;
my $enough_long_pairs=0;
my $num_inv_dups = 0;


open(FH1,"<$first_end") or die;
open(FH2,"<$second_end") or die;

while(my $seq1 = <FH1>, my $seq2 = <FH2>){
	chomp($seq1);
	chomp($seq2);
		if ($linecount_within == 1)
        {
            # Double check the format
            if($seq1 =~ /^\@/ and $seq2 =~ /^\@/){
            	$seq1 =~ s/\/1$/\/2/;
            	if($seq1 ne $seq2){
            		print STDERR "Incorrect FASTQ file \nLine ${linecount_global}: reads $seq1 are not paired and the same order in both of ends fastq files\n";
            	}
            	else{
            		
            	}
            }
            else{
            	print STDERR "Incorrect FASTQ file \nLine ${linecount_global}: $seq1\nMod4 lines should start with \@\n";

            }
  
        }
        elsif ($linecount_within == 2)
        {
			if(length($seq1) >= $first_N and length($seq2) >= $first_N){
					if(&match($seq1, $seq2)){
						$num_inv_dups++;
					}
					$enough_long_pairs++;
			}

        }
        elsif ($linecount_within == 4)
        {
        	$linecount_within = 0;
        }
	$linecount_within++;
	$linecount_global++;
}
close(FH1);
close(FH2);


#my $percentage = sprintf("%.2f",100*$num_inv_dups/$enough_long_pairs);
#print "There are $enough_long_pairs enough long (>=$first_N) pair of reads in the fastq file\n";
#print "There are $num_inv_dups inverted dups reads (${percentage}%) in the fastq file\n";
#$enough_long_pairs *= 2;
## take care!! Here is not the real mapped reads!! Here is the reads that in the fastq file. This is just for USCEC pipeline convenient. 
$linecount_global--;
$linecount_global /=2 ;
my $percentage = sprintf("%.2f",100*$num_inv_dups/($linecount_global/2));
open(OUT,">$output_log") or die;
print OUT "mapped reads=$linecount_global\n";
print OUT "Inverted Read Pairs=$num_inv_dups\n";
print OUT "inverted Pair Percentage=$percentage\n";
close(OUT);

sub match{
	my $seq= shift @_;
	my $query = shift @_;
	my @seqs=split "",$seq;
	my @queries = split "",$query;
	my $mismatch_count = 0;
	for(my $i=0;$i<$first_N;$i++){
		if(&bisulfite_match($seqs[$i],$queries[$i])){
			
		}
		else{
			$mismatch_count++;
		}
		if($mismatch_count > $mismatch){
			return 0;
		}
	}
	return 1;
	
}

sub bisulfite_match{
	my $first = shift @_;
	my $second = shift @_;
	if($first eq $second){
		return 1;
	}
	else{
		if(($first eq 'T' and $second eq 'C') or ($first eq 'G' and $second eq 'A')){
			return 1;
		}
	}
	return 0;
}
