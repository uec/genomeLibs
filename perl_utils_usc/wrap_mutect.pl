#!/usr/bin/perl
use File::Basename;
use Getopt::Long;

GetOptions ("t=s" => \@tumorBams, "n=s" => \@normalBams, "o=s" => \$output, "r=s" => \$ref, "i=s" => \$interval);


my $numcores = `cat /proc/cpuinfo | grep processor -c`;
chomp $numcores;
#$numcores = $numcores / 2 - 1;

my $confidance = 20;
my $bisulfiteRate = 0.9975;
my $minMapQ = 30;
my $minBaseQ = 5;

my $numcores = `cat /proc/cpuinfo | grep processor -c`;
chomp $numcores;
my $ram = 1.3 * $numcores;
$ram = int($ram);
$ram = 20 if $ram > 20;

my $MUTECT = "/home/uec-00/shared/production/software/mutect/default/mutect.jar";
my $JAVA = "/home/uec-00/shared/production/software/java/default_java6/bin/java -Xmx7G";
my $SAMTOOLS = "/home/uec-00/shared/production/software/samtools/samtools";
my $JAVA7 = "/home/uec-00/shared/production/software/java/default_java7/bin/java -Xmx$ram" . "G";
my $GATKSNP = "/home/uec-00/shared/production/software/GATK2/default/GenomeAnalysisTK.jar";


## define required file by provided reference genome
if($ref =~/hg18/){
	$dbsnp="/home/uec-00/shared/production/software/bissnp/genomic_data/dbsnp_135.hg18.sort.vcf";
	$cosmic = "/auto/uec-00/shared/production/genomic-data-misc/COSMIC/hg18_cosmic_v54_120711.vcf";
	$indel_1 = "/home/uec-00/shared/production/software/bissnp/genomic_data/1000G_phase1.indels.hg18.sort.vcf";
	$indel_2 = "/home/uec-00/shared/production/software/bissnp/genomic_data/Mills_and_1000G_gold_standard.indels.hg18.sites.sort.vcf";
}
elsif($ref =~/male.hg19/){
	$dbsnp="/home/uec-00/shared/production/genomic-data-misc/dbsnp/dbsnp_138.hg19.vcf";
	$interval = "/auto/uec-00/shared/production/genomic-data-misc/agilentSureSelect/S04380110_Regions.bed" unless $interval;
	$cosmic = "/auto/uec-00/shared/production/genomic-data-misc/COSMIC/b37_cosmic_v54_120711.fixed.vcf";
	$dbsnp="/home/uec-00/shared/production/genomic-data-misc/dbsnp/dbsnp_138.hg19.vcf";
        $indel_1 = "/auto/uec-00/shared/production/genomic-data-misc/indels/Mills_and_1000G_gold_standard.indels.hg19.vcf";
}
elsif($ref =~/hg19/){
	$indel_1 = "/home/uec-00/shared/production/software/bissnp/genomic_data/1000G_phase1.indels.hg19.sort.vcf";
        $indel_2 = "/home/uec-00/shared/production/software/bissnp/genomic_data/Mills_and_1000G_gold_standard.indels.hg19.sites.sort.vcf";
	$dbsnp = "/home/uec-00/shared/production/genomic-data-misc/dbsnp/dbsnp_135.hg19.sort.vcf";
	$interval = "/auto/uec-00/shared/production/genomic-data-misc/agilentSureSelect/S04380110_Regions.bed" unless $interval;
	$cosmic = "/auto/uec-00/shared/production/genomic-data-misc/COSMIC/b37_cosmic_v54_120711.fixed.vcf";
}
elsif($ref =~/37/ ){
	$dbsnp = "/home/uec-00/shared/production/software/bissnp/genomic_data/dbsnp_135.b37.vcf";
	$interval = "/auto/uec-00/shared/production/genomic-data-misc/agilentSureSelect/S04380110_Regions.bed" unless $interval;
	$cosmic = "/auto/uec-00/shared/production/genomic-data-misc/COSMIC/b37_cosmic_v54_120711.vcf";
        $indel_1 = "/home/uec-00/shared/production/software/bissnp/genomic_data/1000G_phase1.indels.b37.sort.vcf";
        $indel_2 = "/home/uec-00/shared/production/software/bissnp/genomic_data/Mills_and_1000G_gold_standard.indels.b37.sites.sort.vcf";
}
else{
	die "not support this reference genome file yet";
}

push @normalREBams, recalRA(indelRA($_)) for @normalBams;
push @tumorREBams, recalRA(indelRA($_)) for @tumorBams;

&mutect();


sub mutect
{
        createBai($_) for @normalREBams;
        createBai($_) for @tumorREBams;
	runcmd("cp $cosmic .");

        my $cmd .= "$JAVA -jar $MUTECT ";
        $cmd .= "--analysis_type MuTect ";
        $cmd .= "--reference_sequence $ref ";
        $cmd .= "--cosmic " . basename($cosmic) . " ";
        $cmd .= "--dbsnp $dbsnp ";
        $cmd .= "--intervals $interval ";
        $cmd .= "--input_file:normal $_ " for @normalREBams;
        $cmd .= "--input_file:tumor $_ " for @tumorREBams;
        $cmd .= "--out $output\.mutect_stats.txt ";
        $cmd .= "--coverage_file $output\.mutect_coverage.wig";
        runcmd($cmd);
}

sub indelRA
{
        my $inputBam = shift @_;
        my $outputBam = basename($inputBam) . ".realn.bam";
        my $outputIV = basename($inputBam) . ".intervals";
        my $cmd .= "$JAVA7 -jar $GATKSNP -R $ref ";
        $cmd .= "-T RealignerTargetCreator ";
        $cmd .= "-I $inputBam ";
        $cmd .= "-o $outputIV ";
        $cmd .= "-known $indel_1 ";
        $cmd .= "-known $indel_2 " if $indel_2 && $ref !~/mm9/;
        $cmd .= "-nt $numcores ";
        $cmd .= "-U ALLOW_N_CIGAR_READS ";
        runcmd($cmd);

        my $cmd .= "$JAVA7 -jar $GATKSNP -R $ref ";
        $cmd .= "-T IndelRealigner ";
        $cmd .= "-targetIntervals " . basename($inputBam) . ".intervals ";
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
        my $outputBam = basename($inputBam) . ".recal.bam";
        my $outputTable = basename($outputBam) . ".recal_table";
        my $cmd .= "$JAVA7 -jar $GATKSNP -R $ref ";
        $cmd .= "-T BaseRecalibrator ";
        $cmd .= "-I $inputBam ";
        $cmd .= "-o $outputTable ";
        $cmd .= "-knownSites $dbsnp ";
        $cmd .= "-knownSites $indel_1 ";
        $cmd .= "-nct $numcores ";
        $cmd .= "-U ALLOW_N_CIGAR_READS ";
        runcmd($cmd);

        my $cmd .= "$JAVA7 -jar $GATKSNP -R $ref ";
        $cmd .= "-T PrintReads ";
        $cmd .= "-I $inputBam ";
        $cmd .= "-o $outputBam ";
        $cmd .= "-BQSR $outputTable ";
        $cmd .= "-U ALLOW_N_CIGAR_READS ";
        runcmd($cmd);
        return $outputBam;
}

sub createBai
{
	my $bam = shift @_;
	my $bai = $bam;
	$bai =~ s/bam$/bai/;
	runcmd("$SAMTOOLS index $bam") if ! -e "$bam.bai" && ! -e $bai;
}


sub runcmd
{
	my $cmd=shift @_;
	my $caller=(caller(1))[3];
	print STDERR "$caller\t$cmd\n";
	system($cmd);
}