#!/usr/bin/perl
use Cwd 'abs_path';
use File::Basename;
use lib dirname (__FILE__);
use EpigenomeUtils;

my $output = shift @ARGV || die "specify output file to write to ";

my %progs;
$progs{SAMTOOLS} = "$SOFTWAREROOT/software/samtools/samtools";
$progs{PICARD} = "$SOFTWAREROOT/software/picard/default/";
$progs{JAVA} = "$SOFTWAREROOT/software/java/default/bin/java";
$progs{BWA} = "$SOFTWAREROOT/software/bwa/default/";
$progs{TOPHAT} = "$SOFTWAREROOT/software/tophat/default/";
$progs{TOPHAT2} = "$SOFTWAREROOT/software/tophat2/default/";
$progs{BOWTIE} = "$SOFTWAREROOT/software/bowtie/default/";
$progs{BOWTIE2} = "$SOFTWAREROOT/software/bowtie2/default/";
$progs{CUFFLINKS} = "$SOFTWAREROOT/software/cufflinks/default/";
$progs{CUFFLINKS} = "$SOFTWAREROOT/software/cufflinks2/default/";
$progs{UEC_GATK} = "$SOFTWAREROOT/software/uecgatk/default/";
$progs{BSMAP} = "$SOFTWAREROOT/software/bsmap/default/";


open(OUT, ">$output");

for $app (keys %progs)
{
	print OUT "$app=";
	print OUT abs_path($progs{$app});
	print OUT "\n";
}

close OUT;
