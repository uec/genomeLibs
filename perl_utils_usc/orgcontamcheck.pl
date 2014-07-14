#!/usr/bin/perl 
use File::Basename;
# Mon Jul 12 14:03:06 PDT 2010


my $bwa = "/home/uec-00/shared/production/software/perl_utils_usc/wrap_bwa.pl";
my $SAMTOOLS = "/home/uec-00/shared/production/software/samtools/samtools";
my $PICARD = "/home/uec-00/shared/production/software/picard/default";
my $JAVA = "/home/uec-00/shared/production/software/java/default/bin/java";

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



sub runcmd
{
        my $cmd = shift @_;
        print STDERR "$cmd\n";
        system($cmd);
}
