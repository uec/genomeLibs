#!/usr/bin/perl
##This script is used for USC Epigenome center data processing pipeline. Sorted, indexed, have Readgroup tag raw bam file -> BAM with indel realignment -> 
## mark duplicated reads -> BAM with base quality recalibration (Level I) -> bissnp genotyping -> vcf sorting -> vcf filtering, no SB, no SNP cluster (level II) -> 6+2 bed file, strand_combined (level III)

## author: Yaping Liu  lyping1986@gmail.com 
## time: 2012-6-12

#Usege:  bissnp_pipeline.pl bam_file ref_genome_file 

use Getopt::Long;
use File::Basename;

my $skip_indel_align = "";
my $skip_mdups = "";
my $skip_base_recal = "";
my $skip_genotyping = "";
my $skip_vcf_filter = "";
my $skip_vcf_sort = "";
my $skip_filter_bed_file = "";
my $skip_bed_file = "";
my $skip_wig_file = "";

GetOptions( 
	"skip_indel_align" => \$skip_indel_align,
	"skip_mdups" => \$skip_mdups,
	"skip_base_recal" => \$skip_base_recal,
	"skip_genotyping" => \$skip_genotyping,
	"skip_vcf_filter" => \$skip_vcf_filter,
	"skip_vcf_sort" => \$skip_vcf_sort,
	"skip_filter_bed_file" => \$skip_filter_bed_file,
	"skip_bed_file" => \$skip_bed_file,
	"skip_wig_file" => \$skip_wig_file,
);

my $input = shift @ARGV || die "need input bam file";
my $ref = shift @ARGV || die "need reference genome file";

for (@ARGV)
{
	$RRBS = "true" if $_ =~ /RRBS/i;
	$nome_seq_mode = "true" if $_ =~ /NOME/i;
}


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
my $BISSNP = "/home/uec-00/shared/production/software/bissnp/bissnp-default.jar";
my $JAVA = "/home/uec-00/shared/production/software/java/default/bin/java -Xmx$ram" . "G";
my $IGVTOOLS = "/home/uec-00/shared/production/software/igvtools/default/igvtools toTDF";
my $WIG2BW = "/home/uec-00/shared/production/software/UCSC_Browser_Tools/default/wigToBigWig";
my $CHROMSIZE="/home/uec-00/shared/production/software/UCSC_Browser_Tools/default/hg19.chrom.sizes";

my $VCFTOOLS = "/home/uec-00/shared/production/software/bissnp/sortByRefAndCor.pl";
my $VCF2BED = "/home/uec-00/shared/production/software/bissnp/vcf2bed6plus2.pl";
my $VCF2WIG = "/home/uec-00/shared/production/software/bissnp/vcf2wig.pl";
my $VCF2COV = "/home/uec-00/shared/production/software/bissnp/vcf2wig_ct_coverage.pl";
my $VCF2RAWCOV = "/home/uec-00/shared/production/software/bissnp/vcf2wig_raw_coverage.pl";



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
}elsif($ref =~/hg19/){
	$dbsnp="/home/uec-00/shared/production/software/bissnp/genomic_data/dbsnp_135.hg19.sort.vcf";
	$indel_1 = "/home/uec-00/shared/production/software/bissnp/genomic_data/1000G_phase1.indels.hg19.sort.vcf";
	$indel_2 = "/home/uec-00/shared/production/software/bissnp/genomic_data/Mills_and_1000G_gold_standard.indels.hg19.sites.sort.vcf";
	#$ref="/home/uec-00/shared/production/genomes/hg19_rCRSchrm/hg19_rCRSchrm.fa";
	$interval = "/home/uec-00/shared/publicationData/bissnp2011/whole_genome_interval_list.hg19.bed";
	
}elsif($ref =~/37/){
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
}elsif($ref =~/mm10/){
	$dbsnp="/home/uec-00/shared/production/software/bissnp/genomic_data/mgp.v3.snps.rsIDdbSNPv137.mm10.vcf";
	$indel_1 = "/home/uec-00/shared/production/software/bissnp/genomic_data/mgp.v3.indels.rsIDdbSNPv137.mm10.vcf";
	#$ref="/home/uec-00/shared/production/genomes/mm10/mm10.fa";
	$interval = "/home/uec-00/shared/publicationData/bissnp2011/whole_genome_interval_list.mm10.bed";
	$CHROMSIZE="/home/uec-00/shared/production/software/UCSC_Browser_Tools/default/mm10.chrom.sizes";
}
else{
	die "not support this reference genome file yet";
}

## intermdediate file name 
my $input_bam_realign = $input;
$input_bam_realign =~ s/\.bam//;
$input_bam_realign .= ".realign.bam";

my $input_bam_realign_mdups = $input_bam_realign;
$input_bam_realign_mdups =~ s/\.bam//;
$input_bam_realign_mdups .= ".mdups.bam" unless $RRBS;
$input_bam_realign_mdups .= ".withdups.bam" if $RRBS;


my $input_bam_realign_mdups_recal = $input_bam_realign_mdups;
$input_bam_realign_mdups_recal =~ s/\.bam//;
$input_bam_realign_mdups_recal .= ".recal.bam";

my $indel_target_interval = $input;
$indel_target_interval =~ s/\.bam//;
$indel_target_interval .= ".indels.intervals";
my $recalFile_before=$input_bam_realign_mdups_recal;
$recalFile_before =~ s/\.bam//;
$recalFile_before .= ".beforeRecal.txt";
my $recalFile_after=$input_bam_realign_mdups_recal;
$recalFile_after =~ s/\.bam//;
$recalFile_after .= ".afterRecal.txt";

my $tmp_dir="./";
my $metrics_file=$input_bam_realign_mdups.".dupsMetrics.txt";

my $vcf_unsorted_cpg = $input_bam_realign_mdups_recal;
my $vcf_unsorted_snp = $input_bam_realign_mdups_recal;
$vcf_unsorted_cpg =~ s/\.bam//;
if($nome_seq_mode ne ""){
	$vcf_unsorted_cpg .= ".cytosine.raw.vcf";
}
else{
	$vcf_unsorted_cpg .= ".cpg.raw.vcf";
}

$vcf_unsorted_snp =~ s/\.bam//;
$vcf_unsorted_snp .= ".snp.raw.vcf";

my $vcf_sorted_cpg = $vcf_unsorted_cpg;
my $vcf_sorted_snp = $vcf_unsorted_snp;
$vcf_sorted_cpg =~ s/\.vcf//;
$vcf_sorted_cpg .= ".sort.vcf";
$vcf_sorted_snp =~ s/\.vcf//;
$vcf_sorted_snp .= ".sort.vcf";

my $vcf_sorted_cpg_filtered = $vcf_sorted_cpg;
my $vcf_sorted_snp_filtered = $vcf_sorted_snp;
$vcf_sorted_cpg_filtered =~ s/\.raw\.sort\.vcf//;
$vcf_sorted_cpg_filtered .= ".filtered.sort.vcf";
$vcf_sorted_snp_filtered =~ s/\.raw\.sort\.vcf//;
$vcf_sorted_snp_filtered .= ".filtered.sort.vcf";

my $wig_cpg = $vcf_sorted_cpg_filtered;
my $wig_cpg_cov = $vcf_sorted_cpg_filtered;
my $wig_hcg = $vcf_sorted_cpg_filtered;
my $wig_hcg_cov = $vcf_sorted_cpg_filtered;
my $wig_hcg_tdf;
my $wig_hcg_tdf_cov;



$wig_cpg =~ s/\.vcf//;

if($nome_seq_mode ne ""){
	$wig_cpg .= ".GCH.wig";
	$wig_cpg_cov =~ s/.vcf//;
	$wig_cpg_cov .= ".GCH.ct_coverage.wig";
	
	$wig_hcg =~ s/.vcf//;
	$wig_hcg .= ".HCG.wig";
	$wig_hcg_cov =~ s/.vcf//;
	$wig_hcg_cov .= ".HCG.ct_coverage.wig";
	
	$wig_hcg_tdf = $wig_hcg;
	$wig_hcg_tdf =~ s/.wig//;
	$wig_hcg_tdf_cov = $wig_hcg_cov;
	$wig_hcg_tdf_cov =~ s/.wig//;
	$wig_hcg_tdf .= ".bw";
	$wig_hcg_tdf_cov .= ".bw";

}
else{
	$wig_cpg .= ".CG.wig";
	$wig_cpg_cov =~ s/.vcf//;
	$wig_cpg_cov .= ".CG.ct_coverage.wig";
}

my $wig_cpg_raw_cov = $vcf_sorted_cpg;
$wig_cpg_raw_cov =~ s/\.vcf$/.CG.raw_coverage.wig/;
my $wig_tdf_raw_cov = $wig_cpg_raw_cov;
$wig_tdf_raw_cov =~ s/\.wig$/.bw/;

my $wig_tdf = $wig_cpg;
$wig_tdf =~ s/\.wig//;
my $wig_tdf_cov = $wig_cpg_cov;
$wig_tdf_cov =~ s/\.wig//;
my $genome_version;
if($ref =~/hg18/){
	$genome_version = "hg18";
}
elsif($ref =~/hg19/){
	$genome_version = "hg19";
}
elsif($ref =~/37/){
	$genome_version = "b37";
}
elsif($ref =~/mm9/){
	$genome_version = "mm9";
}
elsif($ref =~/mm10/){
	$genome_version = "mm10";
}
$wig_tdf .= ".bw";
$wig_tdf_cov .= ".bw";

unless($skip_indel_align ne "" || -e $input_bam_realign || -e $input_bam_realign_mdups || -e $input_bam_realign_mdups_recal){
	&bam_indel_realign();
}
unless($skip_mdups ne "" || -e $input_bam_realign_mdups || -e $input_bam_realign_mdups_recal){
	&bam_mdups();
}
unless($skip_base_recal ne "" || -e $input_bam_realign_mdups_recal){
	&bam_base_recalibration();
}
my $cmd="rm $input_bam_realign";
system($cmd);
$input_bam_realign=~s/\.bam/.bai/;
$cmd="rm $input_bam_realign";
system($cmd);
$cmd="rm $input_bam_realign_mdups";
system($cmd);
$input_bam_realign_mdups=~s/\.bam/.bai/;
$cmd="rm $input_bam_realign_mdups";
system($cmd);
if($skip_genotyping eq ""){
	&bissnp();
}
if($skip_vcf_sort eq ""){
	&vcf_sort();
}
$cmd="rm $vcf_unsorted_cpg";
system($cmd);
$cmd="rm $vcf_unsorted_snp";
system($cmd);
if($skip_vcf_filter eq ""){
	&vcf_filter();
}
if($skip_bed_file eq "" && $nome_seq_mode eq ""){
	&vcf2bed6plus2();
}
if($skip_filter_bed_file eq ""){
	&vcf2bed6plus2_filter();
}

if($skip_wig_file eq ""){
	&vcf2tdf();
}

sub bam_indel_realign{

	my $header .= "$JAVA -jar $BISSNP -R $ref ";	
	$header .= "-I $input ";
	#$header .= "-T BisulfiteRealignerTargetCreator -L $interval ";
	$header .= "-T BisulfiteRealignerTargetCreator ";
	$header .= "-known $indel_1 ";
	if($ref !~/mm\d+/){
		$header .= "-known $indel_2 ";
	}
	$header .= "-o $indel_target_interval -nt $numcores\n ";
										
	#realign 
	$header .= "$JAVA -jar $BISSNP -R $ref ";	
	$header .= "-I $input ";
	$header .= "-T BisulfiteIndelRealigner -targetIntervals $indel_target_interval ";
	$header .= "-known $indel_1 ";
	if($ref !~/mm\d+/){
	$header .= "-known $indel_2 ";
	}
	$header .= "-compress 5 -cigar ";
	$header .= "-o $input_bam_realign\n ";
	
	system($header)== 0 || die "system call to bam_indel_realign() failed: $?\n" ;
	print "$header\n";
}

#mark dups if not RRBS / unmark if rrbs
sub bam_mdups{
	my $header .= "$JAVA -jar $PICARD/MarkDuplicates.jar I=$input_bam_realign O=$input_bam_realign_mdups ";
	$header .= "CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT MAX_RECORDS_IN_RAM=3000000 ";
	$header .= "METRICS_FILE=$metrics_file ";
	$header .= "TMP_DIR=$tmp_dir \n";
	#handle RRBS
	$header = "$JAVA -jar $PICARD/RevertSam.jar INPUT=$input_bam_realign OUTPUT=$input_bam_realign_mdups VALIDATION_STRINGENCY=SILENT MAX_RECORDS_IN_RAM=3000000 SORT_ORDER=coordinate REMOVE_DUPLICATE_INFORMATION=true REMOVE_ALIGNMENT_INFORMATION=false CREATE_INDEX=true RESTORE_ORIGINAL_QUALITIES=false ATTRIBUTE_TO_CLEAR=null TMP_DIR=$tmp_dir " if $RRBS;
	system($header)== 0 || die "system call to bam_mdups() failed: $?\n" ;;
	print "$header\n";
}

sub bam_base_recalibration{
	##1 countCovariant
	my $header .= "$JAVA -jar $BISSNP -R $ref ";	
	$header .= "-I $input_bam_realign_mdups ";
	$header .= "-T BisulfiteCountCovariates -nt $numcores ";
	$header .= "-knownSites $dbsnp ";
	$header .= "-cov ReadGroupCovariate -cov QualityScoreCovariate -cov CycleCovariate ";
	$header .= "-recalFile $recalFile_before \n ";

##2 TableRecalibration
	$header .= "$JAVA -jar $BISSNP -R $ref ";	
	$header .= "-I $input_bam_realign_mdups ";
	$header .= "-T BisulfiteTableRecalibration ";
	$header .= "-o $input_bam_realign_mdups_recal ";
	$header .= "-recalFile $recalFile_before -maxQ 40 \n ";

##3 countCovariantAfterRecalibrate
	$header .= "$JAVA -jar $BISSNP -R $ref ";	
	$header .= "-I $input_bam_realign_mdups_recal ";
	$header .= "-T BisulfiteCountCovariates -nt $numcores ";
	$header .= "-knownSites $dbsnp ";
	$header .= "-cov ReadGroupCovariate -cov QualityScoreCovariate -cov CycleCovariate ";
	$header .= "-recalFile $recalFile_after \n ";
	
	system($header)== 0 || die "system call to bam_base_recalibration() failed: $?\n" ;;
	print "$header\n";
}

sub bissnp{
	my $header .= "$JAVA -jar $BISSNP -R $ref ";
	$header .= "-I $input_bam_realign_mdups_recal ";									
	$header .= "-D $dbsnp -T BisulfiteGenotyper -vfn1 $vcf_unsorted_cpg -vfn2 $vcf_unsorted_snp ";
	#$header .= "-stand_call_conf $confidance -stand_emit_conf 0 -dt NONE -L $interval -bsRate $bisulfiteRate -loc -1 -nt $numcores ";
	if($nome_seq_mode ne ""){
		$header .= "-out_modes EMIT_VARIANT_AND_CYTOSINES -sm GM ";
	}
	$header .= "-toCoverage 99999 -minConv 0 " if $RRBS;
	print STDERR "using bissnp RBBS options\n " if $RRBS;
	$header .= "-stand_call_conf $confidance -stand_emit_conf 0 -dt NONE -bsRate $bisulfiteRate -loc -1 -nt $numcores ";
	$header .= "-minConv 1 -vcfCache 1000000 ";
	$header .= "-mmq $minMapQ ";
	$header .= "-mbq $minBaseQ\n";
	
	system($header)== 0 || die "system call to bissnp() failed: $?\n" ;;
	print "$header\n";
}
 
sub vcf_sort{
	my $header .= "perl $VCFTOOLS --k 1 --c 2 ";
	$header .= "--tmp ".dirname($vcf_unsorted_cpg);
	$header .= " $vcf_unsorted_cpg ";
	$header .= "$ref";
	$header .= ".fai ";
	$header .= "> $vcf_sorted_cpg \n";
	$header .= "perl $VCFTOOLS --k 1 --c 2 ";
	$header .= "--tmp ".dirname($vcf_unsorted_snp)." ";
	$header .= "$vcf_unsorted_snp ";
	$header .= "$ref";
	$header .= ".fai ";
	$header .= "> $vcf_sorted_snp \n";
	
	system($header)== 0 || die "system call to vcf_sort() failed: $?\n" ;;
	print "$header\n";
}

sub vcf_filter{
	my $header .= "$JAVA -jar $BISSNP -R $ref -T VCFpostprocess ";
	#$header .= "-qual $confidance -L $interval ";
	if($nome_seq_mode ne ""){
		$header .= "-C GCH -C HCH -C GCG -C HCG ";
	}
	$header .= "-qual $confidance ";
	$header .= "-oldVcf $vcf_sorted_cpg ";
	$header .= "-snpVcf $vcf_sorted_snp ";
	$header .= "-newVcf $vcf_sorted_cpg_filtered ";
	$header .= "-o $vcf_sorted_cpg_filtered.cpgSummary.txt \n";
	
	$header .= "$JAVA -jar $BISSNP -R $ref -T VCFpostprocess ";
	#$header .= "-qual $confidance -L $interval ";
	if($nome_seq_mode ne ""){
		$header .= "-C GCH -C HCH -C GCG -C HCG ";
	}
	$header .= "-qual $confidance ";
	$header .= "-oldVcf $vcf_sorted_snp ";
	$header .= "-snpVcf $vcf_sorted_snp ";
	$header .= "-newVcf $vcf_sorted_snp_filtered ";
	$header .= "-o $vcf_sorted_snp_filtered.cpgSummary.txt \n";
	
	system($header)== 0 || die "system call to vcf_filter() failed: $?\n" ;;
	print "$header\n";
	
}

##output all CGs, but in seperated strand
sub vcf2bed6plus2{
	my $header .= "perl $VCF2BED --seperate_strand $vcf_sorted_cpg ";
	if($nome_seq_mode ne ""){
		$header .= "GCH \n";
		system($header)== 0 || die "system call to vcf2bed6plus2() failed: $?\n" ;;
		print "$header\n";
		$header = "perl $VCF2BED --seperate_strand $vcf_sorted_cpg HCG \n";
	}
	else{
		$header .= "CG \n";
	}
	
	system($header)== 0 || die "system call to vcf2bed6plus2() failed: $?\n" ;;
	print "$header\n";
}

##output only good CGs{, for NOMe-seq, it will not combine both of the strand..
sub vcf2bed6plus2_filter{
	my $header .= "perl $VCF2BED --only_good_call $vcf_sorted_cpg_filtered ";
	if($nome_seq_mode ne ""){
		$header .= "GCH --seperate_strand \n";
		system($header)== 0 || die "system call to vcf2bed6plus2_filter() failed: $?\n" ;;
		print "$header\n";
		$header = "perl $VCF2BED --seperate_strand --only_good_call $vcf_sorted_cpg_filtered HCG \n";
	}
	else{
		$header .= "CG \n";
	}
	
	system($header)== 0 || die "system call to vcf2bed6plus2_filter() failed: $?\n" ;;
	print "$header\n";
}

sub vcf2tdf{
	my $header = "perl $VCF2WIG $vcf_sorted_cpg_filtered ";
	if($nome_seq_mode ne ""){
		$header .= "GCH \n";
		runcmd($header,"vcf2tdf()");
		$header = "perl $VCF2WIG $vcf_sorted_cpg_filtered HCG \n";
	}
	else{
		$header .= "CG \n";
	}
	runcmd($header,"vcf2tdf()");
	
	$header = "perl $VCF2COV $vcf_sorted_cpg_filtered ";
	if($nome_seq_mode ne ""){
		$header .= "GCH \n";
		runcmd($header,"vcf2tdf()");
		$header = "perl $VCF2COV $vcf_sorted_cpg_filtered HCG \n";
	}
	else{
		$header .= "CG \n";
	}
	runcmd($header,"vcf2tdf()");
	
	if($nome_seq_mode eq ""){	
		$header = "perl $VCF2RAWCOV $vcf_sorted_cpg CG \n";
		runcmd($header,"vcf2tdf()");
	}
	#$header = "$IGVTOOLS $wig_cpg $wig_tdf $genome_version \n";
	$header = "$WIG2BW $wig_cpg $CHROMSIZE $wig_tdf\n";
	runcmd($header,"vcf2tdf()");
	#$header = "$IGVTOOLS $wig_cpg_cov $wig_tdf_cov $genome_version \n";	
	$header = "$WIG2BW $wig_cpg_cov $CHROMSIZE $wig_tdf_cov \n";	
	
	runcmd($header,"vcf2tdf()");
	#$header = "rm $wig_cpg\n";
	#$header .= "rm $wig_cpg_cov\n";
	#runcmd($header,"vcf2tdf()");
	
	if($nome_seq_mode eq ""){
		#$header = "$IGVTOOLS $wig_cpg_raw_cov $wig_tdf_raw_cov $genome_version \n";
		$header = "$WIG2BW $wig_cpg_raw_cov $CHROMSIZE $wig_tdf_raw_cov \n";
		runcmd($header,"vcf2tdf()");
		#$header = "rm $wig_cpg_raw_cov\n";
		#runcmd($header,"vcf2tdf()");
	}
	
	runcmd($header,"vcf2tdf()");
	if($nome_seq_mode ne ""){
		#$header = "$IGVTOOLS $wig_hcg $wig_hcg_tdf $genome_version \n";
		$header = "$WIG2BW $wig_hcg $CHROMSIZE $wig_hcg_tdf \n";
		runcmd($header,"vcf2tdf()");
		#$header = "$IGVTOOLS $wig_hcg_cov $wig_hcg_tdf_cov $genome_version \n";
		$header = "$WIG2BW $wig_hcg_cov $CHROMSIZE $wig_hcg_tdf_cov \n";
		
		runcmd($header,"vcf2tdf()");
		
		#$header = "rm $wig_hcg\n";
		#$header = "rm $wig_hcg_cov\n";
		#runcmd($header,"vcf2tdf()");
		
	}
	
}

sub runcmd{
	my $cmd=shift @_;
	my $module_name=shift @_;
	system($cmd)== 0 || die "system call to $module_name failed: $?\n" ;
	print "$cmd\n";
}
