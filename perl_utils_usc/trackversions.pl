#!/usr/bin/perl
use Cwd 'abs_path';

my $output = shift @ARGV || die "specify output file to write to ";

my %progs;
$progs{SAMTOOLS} = "/home/uec-00/shared/production/software/samtools/samtools";
$progs{PICARD} = "/home/uec-00/shared/production/software/picard/default/";
$progs{JAVA} = "/home/uec-00/shared/production/software/java/default/bin/java";
$progs{BWA} = "/home/uec-00/shared/production/software/bwa/default/";
$progs{TOPHAT} = "/home/uec-00/shared/production/software/tophat/default/";
$progs{BOWTIE} = "/home/uec-00/shared/production/software/bowtie/default/";
$progs{CUFFLINKS} = "/home/uec-00/shared/production/software/cufflinks/default/";
$progs{UEC_GATK} = "/home/uec-00/shared/production/software/uecgatk/default/";
$progs{BSMAP} = "/home/uec-00/shared/production/software/bsmap/default/";


open(OUT, ">$output");

for $app (keys %progs)
{
        print OUT "$app=";
        print OUT abs_path($progs{$app});
        print OUT "\n";
}

close OUT;
