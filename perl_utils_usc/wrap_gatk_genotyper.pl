#!/usr/bin/perl
use File::Basename;
use lib dirname (__FILE__);
use EpigenomeUtils;

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
$ram = int($ram);
$ram = 20 if $ram > 20;


my $IGVTOOLS = "$SOFTWAREROOT/igvtools/default/igvtools toTDF";
my $WIG2BW = "$SOFTWAREROOT/UCSC_Browser_Tools/default/wigToBigWig";
my $CHROMSIZE="$SOFTWAREROOT/UCSC_Browser_Tools/default/hg19.chrom.sizes";

my $VCFTOOLS = "$SOFTWAREROOT/bissnp/sortByRefAndCor.pl";
my $VCF2BED = "$SOFTWAREROOT/bissnp/vcf2bed6plus2.pl";
my $VCF2WIG = "$SOFTWAREROOT/bissnp/vcf2wig.pl";
my $VCF2COV = "$SOFTWAREROOT/bissnp/vcf2wig_ct_coverage.pl";
my $VCF2RAWCOV = "$SOFTWAREROOT/bissnp/vcf2wig_raw_coverage.pl";

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
	$dbsnp="$SOFTWAREROOT/bissnp/genomic_data/dbsnp_135.hg18.sort.vcf";
	$indel_1 = "$SOFTWAREROOT/bissnp/genomic_data/1000G_phase1.indels.hg18.sort.vcf";
	$indel_2 = "$SOFTWAREROOT/bissnp/genomic_data/Mills_and_1000G_gold_standard.indels.hg18.sites.sort.vcf";
	#$ref="/home/uec-00/shared/production/genomes/hg18_unmasked/hg18_unmasked.plusContam.fa";
	$interval = "$PUBLICATIONDATA/bissnp2011/whole_genome_interval_list.hg18.bed";
	$CHROMSIZE="$SOFTWAREROOT/UCSC_Browser_Tools/default/hg18.chrom.sizes";
}
elsif($ref =~/male.hg19/){
	$dbsnp="$GENOMEROOT/genomic-data-misc/dbsnp/dbsnp_138.hg19.vcf";
	$indel_1 = "$GENOMEROOT/genomic-data-misc/indels/Mills_and_1000G_gold_standard.indels.hg19.vcf";
	#$indel_1 = "$SOFTWAREROOT/bissnp/genomic_data/1000G_phase1.indels.hg19.sort.vcf";
	#$indel_2 = "$SOFTWAREROOT/bissnp/genomic_data/Mills_and_1000G_gold_standard.indels.hg19.sites.sort.vcf";
	$interval = "$PUBLICATIONDATA/bissnp2011/whole_genome_interval_list.hg19.bed";
}	
elsif($ref =~/hg19/){
	$dbsnp="$SOFTWAREROOT/bissnp/genomic_data/dbsnp_135.hg19.sort.vcf";
	$indel_1 = "$SOFTWAREROOT/bissnp/genomic_data/1000G_phase1.indels.hg19.sort.vcf";
	$indel_2 = "$SOFTWAREROOT/bissnp/genomic_data/Mills_and_1000G_gold_standard.indels.hg19.sites.sort.vcf";
	#$ref="/home/uec-00/shared/production/genomes/hg19_rCRSchrm/hg19_rCRSchrm.fa";
	$interval = "$PUBLICATIONDATA/bissnp2011/whole_genome_interval_list.hg19.bed";
}	
elsif($ref =~/37/ ){
	$dbsnp="$SOFTWAREROOT/bissnp/genomic_data/dbsnp_135.b37.vcf";
	$indel_1 = "$SOFTWAREROOT/bissnp/genomic_data/1000G_phase1.indels.b37.sort.vcf";
	$indel_2 = "$SOFTWAREROOT/bissnp/genomic_data/Mills_and_1000G_gold_standard.indels.b37.sites.sort.vcf";
	#$ref="$SOFTWAREROOT/bissnp/genomic_data/GRCh37-lite.fa";
	$interval = "$PUBLICATIONDATA/bissnp2011/whole_genome_interval_list.hg19.bed";
	$CHROMSIZE="/export/uec-gs1/laird/users/yaping/data/genome_data/genome_interval/GRCh37.chrom.sizes";

}elsif($ref =~/mm9/){
	$dbsnp="$SOFTWAREROOT/bissnp/genomic_data/mouse-20111102-snps-all.annotated.mm9.vcf";
	$indel_1 = "$SOFTWAREROOT/bissnp/genomic_data/mouse-20110602-callable-dinox-indels.annot.mm9.vcf";
	#$ref="/home/uec-00/shared/production/genomes/mm9_unmasked/mm9_unmasked.fa";
	$interval = "$PUBLICATIONDATA/bissnp2011/whole_genome_interval_list.mm9.bed";
	$CHROMSIZE="$SOFTWAREROOT/UCSC_Browser_Tools/default/mm9.chrom.sizes";
}
else{
	die "not support this reference genome file yet";
}

my $indelBam = indelRA($input);
my $recalBam = recalRA($indelBam);
my $vcf = gatksnp($recalBam);
my $vcfMetrics = gatksnpVariantEval($vcf);


sub indelRA
{
        my $inputBam = shift @_;
        my $outputBam = $inputBam . ".realn.bam";
        my $outputIV = $inputBam . ".intervals";
        my $cmd .= "$JAVA -jar $GATKSNP -R $ref ";
        $cmd .= "-T RealignerTargetCreator ";
        $cmd .= "-I $inputBam ";
        $cmd .= "-o $outputIV ";
        $cmd .= "-known $indel_1 ";
        $cmd .= "-known $indel_2 " if $indel_2 && $ref !~/mm9/;
        $cmd .= "-nt $numcores ";
        $cmd .= "-U ALLOW_N_CIGAR_READS ";
        runcmd($cmd);

        my $cmd .= "$JAVA -jar $GATKSNP -R $ref ";
        $cmd .= "-T IndelRealigner ";
        $cmd .= "-targetIntervals $inputBam\.intervals ";
        $cmd .= "-I $inputBam ";
        $cmd .= "-o $outputBam ";
        $cmd .= "-compress 0 ";
        $cmd .= "-known $indel_1 ";
        $cmd .= "-known $indel_2 " if $indel_2 && $ref !~/mm9/;
        $cmd .= "-U ALLOW_N_CIGAR_READS ";
        runcmd($cmd);
        return $outputBam;
}

sub recalRA
{
        my $inputBam = shift @_;
        my $outputBam = $inputBam . ".recal.bam";
        my $outputTable = $outputBam . ".recal_table";
        my $cmd .= "$JAVA -jar $GATKSNP -R $ref ";
        $cmd .= "-T BaseRecalibrator ";
        $cmd .= "-I $inputBam ";
        $cmd .= "-o $outputTable ";
        $cmd .= "-knownSites $dbsnp ";
        $cmd .= "-knownSites $indel_1 ";
        $cmd .= "-nct $numcores ";
        $cmd .= "-U ALLOW_N_CIGAR_READS ";
        runcmd($cmd);

        my $cmd .= "$JAVA -jar $GATKSNP -R $ref ";
        $cmd .= "-T PrintReads ";
        $cmd .= "-I $inputBam ";
        $cmd .= "-o $outputBam ";
        $cmd .= "-BQSR $outputTable ";
        $cmd .= "-U ALLOW_N_CIGAR_READS ";
        runcmd($cmd);
        return $outputBam;
}



sub gatksnp{
	my $inputBam = shift @_;
	my $outputVCF = $input . ".snps.raw.vcf";
	my $cmd .= "$JAVA -jar $GATKSNP -R $ref ";
	$cmd .= "-I $inputBam ";									
	$cmd .= "-T HaplotypeCaller ";
	$cmd .= "--dbsnp $dbsnp ";
	$cmd .= "-o $outputVCF -U ALLOW_N_CIGAR_READS ";
	$cmd .= "-stand_call_conf $confidance -stand_emit_conf 0 -nct $numcores ";
	#$cmd .= "-mmq $minMapQ -mbq $minBaseQ -out_mode EMIT_VARIANTS_ONLY\n";
	#$cmd .= "-mbq $minBaseQ -out_mode EMIT_VARIANTS_ONLY\n";
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
