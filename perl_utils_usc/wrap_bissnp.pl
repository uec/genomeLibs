#!/usr/bin/perl
use strict;

my $numcores = `cat /proc/cpuinfo | grep processor -c`;
chomp $numcores;
#$numcores = $numcores / 2 - 1;

#my $ram = 4 * $numcores;
my $ram = 12;

my $input = $ARGV[0] || die "need input file";
my $bissnp = "/home/uec-00/shared/production/software/bissnp/bissnp-default.jar";
my $vcf = "/home/uec-00/shared/production/software/bissnp/dbsnp_135.hg19.sort.vcf";
my $ref = "/home/uec-00/shared/production/genomes/hg19_rCRSchrm/hg19_rCRSchrm.fa";
my $interval = "/home/uec-00/shared/production/software/bissnp/wholegenome_interval_list.hg19.bed";
my $JAVA = "/home/uec-00/shared/production/software/java/default/bin/java -Xmx$ram" . "G";



my $cmd = "$JAVA -jar $bissnp -R $ref -T BisulfiteGenotyper -I $input -D $vcf -vfn1 $input.cpg.raw.vcf -vfn2 $input.snp.raw.vcf -stand_call_conf 20 -stand_emit_conf 0 -L $interval -out_modes DEFAULT_FOR_TCGA -nt $numcores -rgv hg19 -mbq 0 -mmq 30";

runcmd("$cmd");

sub runcmd
{
        my $cmd = shift @_;
        print STDERR "$cmd\n";
        system($cmd);
}
