#!/usr/bin/perl 
use File::Basename;
use lib dirname (__FILE__);
use EpigenomeUtils;

my $bwa = "$SOFTWAREROOT/perl_utils_usc/wrap_bwa.pl";
my $input = shift @ARGV;
my $sampleSize = shift @ARGV;
$lines = $sampleSize * 4;

my $inputBase = basename($input);

runcmd("head -n $lines $input > $inputBase\.sampled.txt") unless ($input =~ /gz$/);
runcmd("zcat $input | head -n $lines > $inputBase\.sampled.txt") if($input =~ /gz$/);

for my $genome ( @ARGV )
{
	my $genomeBase = basename($genome);
	
	runcmd("$bwa $genome $inputBase\.sampled.txt $inputBase\.$genomeBase\.bam");
	runcmd("$SAMTOOLS flagstat  $inputBase\.$genomeBase\.bam > ContamCheck\.$genomeBase\.$sampleSize\.$inputBase\.flagstat.txt");
}
