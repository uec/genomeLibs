#!/usr/bin/perl
use File::Basename;

my $input = $ARGV[0] || die "need input bam file";
my $ref = $ARGV[1] || die "need reference genome file";

my $numcores = `cat /proc/cpuinfo | grep processor -c`;
chomp $numcores;
#$numcores = $numcores / 2 - 1;

my $confidance = 20;
my $bisulfiteRate = 0.9975;
my $minMapQ = 30;
my $minBaseQ = 5;

my $ram = 1.3 * $numcores;
my $ram = int($ram);
my $SAMTOOLS = "/home/uec-00/shared/production/software/samtools/samtools";
my $PICARD = "/home/uec-00/shared/production/software/picard/default/";
my $GATKSNP = "/home/uec-00/shared/production/software/GATK2/default/GenomeAnalysisTK.jar";
my $JAVA = "/home/uec-00/shared/production/software/java/default_java7/bin/java -Xmx$ram" . "G";
my $IGVTOOLS = "/home/uec-00/shared/production/software/igvtools/default/igvtools toTDF";
my $WIG2BW = "/home/uec-00/shared/production/software/UCSC_Browser_Tools/default/wigToBigWig";
my $CHROMSIZE="/home/uec-00/shared/production/software/UCSC_Browser_Tools/default/hg19.chrom.sizes";

my $VCFTOOLS = "/home/uec-00/shared/production/software/bissnp/sortByRefAndCor.pl";
my $VCF2BED = "/home/uec-00/shared/production/software/bissnp/vcf2bed6plus2.pl";
my $VCF2WIG = "/home/uec-00/shared/production/software/bissnp/vcf2wig.pl";
my $VCF2COV = "/home/uec-00/shared/production/software/bissnp/vcf2wig_ct_coverage.pl";
my $VCF2RAWCOV = "/home/uec-00/shared/production/software/bissnp/vcf2wig_raw_coverage.pl";

my $confidance = 20;
my $bisulfiteRate = 0.9975;
my $minMapQ = 30;
my $minBaseQ = 5;


my $dbsnp;
my $indel_1;
my $indel_2;
my $interval;


## define required file by provided reference genome
if($ref =~/hg18/){
	$dbsnp="/home/uec-00/shared/production/software/bissnp/genomic_data/dbsnp_135.hg18.sort.vcf";
	$indel_1 = "/home/uec-00/shared/production/software/bissnp/genomic_data/1000G_phase1.indels.hg18.sort.vcf";
	$indel_2 = "/home/uec-00/shared/production/software/bissnp/genomic_data/Mills_and_1000G_gold_standard.indels.hg18.sites.sort.vcf";
	#$ref="/home/uec-00/shared/production/genomes/hg18_unmasked/hg18_unmasked.plusContam.fa";
	$interval = "/home/uec-00/shared/publicationData/bissnp2011/whole_genome_interval_list.hg18.bed";
	$CHROMSIZE="/home/uec-00/shared/production/software/UCSC_Browser_Tools/default/hg18.chrom.sizes";
}
elsif($ref =~/male.hg19/){
	$dbsnp="/home/uec-00/shared/production/genomic-data-misc/dbsnp/dbsnp_138.hg19.vcf";
	$indel_1 = "/home/uec-00/shared/production/software/bissnp/genomic_data/1000G_phase1.indels.hg19.sort.vcf";
	$indel_2 = "/home/uec-00/shared/production/software/bissnp/genomic_data/Mills_and_1000G_gold_standard.indels.hg19.sites.sort.vcf";
	$interval = "/home/uec-00/shared/publicationData/bissnp2011/whole_genome_interval_list.hg19.bed";
}	
elsif($ref =~/hg19/){
	$dbsnp="/home/uec-00/shared/production/software/bissnp/genomic_data/dbsnp_135.hg19.sort.vcf";
	$indel_1 = "/home/uec-00/shared/production/software/bissnp/genomic_data/1000G_phase1.indels.hg19.sort.vcf";
	$indel_2 = "/home/uec-00/shared/production/software/bissnp/genomic_data/Mills_and_1000G_gold_standard.indels.hg19.sites.sort.vcf";
	#$ref="/home/uec-00/shared/production/genomes/hg19_rCRSchrm/hg19_rCRSchrm.fa";
	$interval = "/home/uec-00/shared/publicationData/bissnp2011/whole_genome_interval_list.hg19.bed";
}	
elsif($ref =~/37/ ){
	$dbsnp="/home/uec-00/shared/production/software/bissnp/genomic_data/dbsnp_135.b37.vcf";
	$indel_1 = "/home/uec-00/shared/production/software/bissnp/genomic_data/1000G_phase1.indels.b37.sort.vcf";
	$indel_2 = "/home/uec-00/shared/production/software/bissnp/genomic_data/Mills_and_1000G_gold_standard.indels.b37.sites.sort.vcf";
	#$ref="/home/uec-00/shared/production/software/bissnp/genomic_data/GRCh37-lite.fa";
	$interval = "/home/uec-00/shared/publicationData/bissnp2011/whole_genome_interval_list.hg19.bed";
	$CHROMSIZE="/export/uec-gs1/laird/users/yaping/data/genome_data/genome_interval/GRCh37.chrom.sizes";

}elsif($ref =~/mm9/){
	$dbsnp="/home/uec-00/shared/production/software/bissnp/genomic_data/mouse-20111102-snps-all.annotated.mm9.vcf";
	$indel_1 = "/home/uec-00/shared/production/software/bissnp/genomic_data/mouse-20110602-callable-dinox-indels.annot.mm9.vcf";
	#$ref="/home/uec-00/shared/production/genomes/mm9_unmasked/mm9_unmasked.fa";
	$interval = "/home/uec-00/shared/publicationData/bissnp2011/whole_genome_interval_list.mm9.bed";
	$CHROMSIZE="/home/uec-00/shared/production/software/UCSC_Browser_Tools/default/mm9.chrom.sizes";
}
else{
	die "not support this reference genome file yet";
}

my $vcf = gatksnp($input);
my $vcfMetrics = gatksnpVariantEval($vcf);

sub gatksnp{
	my $inputBam = shift @_;
	my $outputVCF = $inputBam . ".snps.raw.vcf";
	my $cmd .= "$JAVA -jar $GATKSNP -R $ref ";
	$cmd .= "-I $inputBam ";									
	$cmd .= "-D $dbsnp -T UnifiedGenotyper -o $outputVCF ";
	$cmd .= "-stand_call_conf $confidance -stand_emit_conf 0 -nt $numcores ";
	#$cmd .= "-mmq $minMapQ -mbq $minBaseQ -out_mode EMIT_VARIANTS_ONLY\n";
	$cmd .= "-mbq $minBaseQ -out_mode EMIT_VARIANTS_ONLY\n";
	runcmd($cmd);
	return $outputVCF;
}

sub gatksnpVariantEval{
	my $inputVCF = shift @_;
	my $outputTxt = $inputVCF . ".varianteval.metric.txt";
	my $cmd .= "$JAVA -jar $GATKSNP -R $ref ";
	$cmd .= "--eval $inputVCF ";									
	$cmd .= "-D $dbsnp -T VariantEval -o $outputTxt ";
	$cmd .= "-nt $numcores";
	runcmd($cmd);
	return $outputTxt;
}

sub vcf_sort{
	my $cmd .= "perl $VCFTOOLS --k 1 --c 2 ";
	$cmd .= "--tmp ".dirname($vcf_unsorted_cpg);
	$cmd .= " $vcf_unsorted_cpg ";
	$cmd .= "$ref";
	$cmd .= ".fai ";
	$cmd .= "> $vcf_sorted_cpg \n";
	$cmd .= "perl $VCFTOOLS --k 1 --c 2 ";
	$cmd .= "--tmp ".dirname($vcf_unsorted_snp)." ";
	$cmd .= "$vcf_unsorted_snp ";
	$cmd .= "$ref";
	$cmd .= ".fai ";
	$cmd .= "> $vcf_sorted_snp \n";
	
	system($cmd)== 0 || die "system call to vcf_sort() failed: $?\n" ;;
	print "$cmd\n";
}

sub runcmd{
	my $cmd=shift @_;
	my $caller=(caller(1))[3];
	print STDERR "$caller\t$cmd\n";
	system($cmd);
}
