#!/usr/bin/perl

## author: Yaping Liu  lyping1986@gmail.com 
## time: 2013-2-19

#Usege:  bissnp_trinuc_sample.pl $outputfile bam_file ref_genome_file 

use strict;
use Getopt::Long;
my $use_bad_mates;

my $metric_unsorted_cpg = $ARGV[0] || die "specify output txt file";
my $input = $ARGV[1] || die "need input bam file";
my $ref = $ARGV[2] || die "need reference genome file";
my $interval;

my $numcores = `cat /proc/cpuinfo | grep processor -c`;
chomp $numcores;
#$numcores = $numcores / 2 - 1;
#my $numcores = 4;
my $confidance = 20;
my $bisulfiteRate = 0.9975;
my $minMapQ = 30;
my $minBaseQ = 5;

#my $ram = 4 * $numcores;
my $ram = 12;

my $BISSNP = "/home/uec-00/shared/production/software/bissnp/bissnp-default.jar";
my $JAVA = "/home/uec-00/shared/production/software/java/default/bin/java -Xmx$ram" . "G";


my $vcf_unsorted_cpg=$input;
$vcf_unsorted_cpg =~ s/\.bam$/.cytosine.raw.vcf/;
my $vcf_summary=$vcf_unsorted_cpg.".MethySummarizeList.txt";
my $dbsnp;
my @characters=("A","C","G","T");



## define required file by provided reference genome
if($ref =~/hg18/){
	$dbsnp="/home/uec-00/shared/production/software/bissnp/genomic_data/dbsnp_135.hg18.sort.vcf";
	$interval = $ARGV[3] || "chr21";

}
elsif($ref =~/hg19/){
	$dbsnp="/home/uec-00/shared/production/software/bissnp/genomic_data/dbsnp_135.hg19.sort.vcf";
	$interval = $ARGV[3] || "chr21";
	
}
elsif($ref =~/NC_001416/){
	undef $interval;
	undef $dbsnp;
	$use_bad_mates = 1;
}
elsif($ref =~/37/){
	$dbsnp="/home/uec-00/shared/production/software/bissnp/genomic_data/dbsnp_135.b37.vcf";
	$interval = $ARGV[3] || "chr21";
}
elsif($ref =~/mm9/){
	$dbsnp="/home/uec-00/shared/production/software/bissnp/genomic_data/mouse-20111102-snps-all.annotated.mm9.vcf";
	$interval = $ARGV[3] || "chr16";
}

&bissnp();
&output_tri_file();
my $cmd ="rm $vcf_unsorted_cpg\n";
$cmd .="rm $vcf_summary\n";
system($cmd);
print "$cmd\n";


sub bissnp{
	my $header .= "$JAVA -jar $BISSNP -R $ref ";
	$header .= "-I $input ";									
	$header .= "-T BisulfiteGenotyper -vfn1 $vcf_unsorted_cpg ";
	
	$header .= "-L $interval " if $interval;
	$header .= "-D $dbsnp " if $dbsnp;
	$header .= "-out_modes EMIT_ALL_CYTOSINES -sm BM ";
	
	my $c_str="";
	foreach my $a(@characters){
		#$c_str.="";
		foreach my $b(@characters){
			$c_str.="-C ${a}C${b},2 ";	
		}
	}
	$header .= "$c_str";
	$header .= "-badMate " if $use_bad_mates;
	$header .= "-stand_call_conf $confidance -stand_emit_conf 0 -dt NONE -bsRate $bisulfiteRate -nt $numcores ";
	$header .= "-minConv 1 -vcfCache 1000000 ";
	$header .= "-mmq $minMapQ ";
	$header .= "-mbq $minBaseQ\n";
	
	system($header);
	print "$header\n";
}


sub output_tri_file{
	open(FH,"<$vcf_summary") or die "can't open $vcf_summary";
	open(OUT,">$metric_unsorted_cpg") or die "can't open $metric_unsorted_cpg";
	my $out_flag=-1;
	while(<FH>){
		chomp;
		print OUT "$_\n" if($out_flag==1);
		$out_flag=1 if($_=~/^##Methylation summary in Read Group:/);
	}
	close(FH);
	close(OUT);
}

