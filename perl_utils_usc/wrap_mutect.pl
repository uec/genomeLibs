#!/usr/bin/perl
use File::Basename;
use Getopt::Long;

GetOptions ("t=s" => \@tumorBams, "n=s" => \@normalBams, "o=s" => \$output, "r=s" => \$ref);


my $numcores = `cat /proc/cpuinfo | grep processor -c`;
chomp $numcores;
#$numcores = $numcores / 2 - 1;

my $confidance = 20;
my $bisulfiteRate = 0.9975;
my $minMapQ = 30;
my $minBaseQ = 5;

my $MUTECT = "/home/uec-00/shared/production/software/mutect/default/mutect.jar";
my $JAVA = "/home/uec-00/shared/production/software/java/default_java6/bin/java -Xmx7G";
my $SAMTOOLS = "/home/uec-00/shared/production/software/samtools/samtools";


## define required file by provided reference genome
if($ref =~/hg18/){
	$dbsnp="/home/uec-00/shared/production/software/bissnp/genomic_data/dbsnp_135.hg18.sort.vcf";
	$cosmic = "/auto/uec-00/shared/production/genomic-data-misc/COSMIC/hg18_cosmic_v54_120711.vcf";
}
elsif($ref =~/male.hg19/){
	$dbsnp="/home/uec-00/shared/production/genomic-data-misc/dbsnp/dbsnp_138.hg19.vcf";
	$interval = "/auto/uec-00/shared/production/genomic-data-misc/agilentSureSelect/S04380110_Regions.bed";
	$cosmic = "/auto/uec-00/shared/production/genomic-data-misc/COSMIC/b37_cosmic_v54_120711.fixed.vcf";
}	
elsif($ref =~/hg19/){
	$dbsnp="/home/uec-00/shared/production/genomic-data-misc/dbsnp/dbsnp_135.hg19.sort.vcf";
	$interval = "/auto/uec-00/shared/production/genomic-data-misc/agilentSureSelect/S04380110_Regions.bed";
	$cosmic = "/auto/uec-00/shared/production/genomic-data-misc/COSMIC/b37_cosmic_v54_120711.fixed.vcf";
}	
elsif($ref =~/37/ ){
	$dbsnp="/home/uec-00/shared/production/software/bissnp/genomic_data/dbsnp_135.b37.vcf";
	#$ref="/home/uec-00/shared/production/software/bissnp/genomic_data/GRCh37-lite.fa";
	$interval = "/auto/uec-00/shared/production/genomic-data-misc/agilentSureSelect/S04380110_Regions.bed";
	$cosmic = "/auto/uec-00/shared/production/genomic-data-misc/COSMIC/b37_cosmic_v54_120711.vcf";
}
else{
	die "not support this reference genome file yet";
}

&mutect();


sub mutect
{
        createBai($_) for @normalBams;
        createBai($_) for @tumorBams;

        my $cmd .= "$JAVA -jar $MUTECT ";
        $cmd .= "--analysis_type MuTect ";
        $cmd .= "--reference_sequence $ref ";
        $cmd .= "--cosmic $cosmic ";
        $cmd .= "--dbsnp $dbsnp ";
        $cmd .= "--intervals $interval ";
        $cmd .= "--input_file:normal $_ " for @normalBams;
        $cmd .= "--input_file:tumor $_ " for @tumorBams;
        $cmd .= "--out $output\.mutect_stats.txt ";
        $cmd .= "--coverage_file $output\.mutect_coverage.wig";
        runcmd($cmd);
}

sub createBai
{
	my $bam = shift @_;
	my $bai = $bam;
	$bai =~ s/bam$/bai/;
	runcmd("$SAMTOOLS index $bam") if ! -e "$bam.bai" && ! -e $bai;
}


sub runcmd{
	my $cmd=shift @_;
	my $caller=(caller(1))[3];
	print STDERR "$caller\t$cmd\n";
	system($cmd);
}
