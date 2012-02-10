#!/usr/bin/perl
use strict;

my $numcores = `cat /proc/cpuinfo | grep processor -c`;
chomp $numcores;
$numcores = $numcores / 2 - 1;

my $ram = 4 * $numcores;

my $input = $ARGV[0] || die;
my $bissnp = "/home/uec-00/shared/production/software/bissnp/bissnp-default.jar";
my $rod = "/home/uec-00/shared/production/software/bissnp/dbsnp_132.hg19.sort.rod";
my $ref = "/home/uec-00/shared/production/genomes/hg19_rCRSchrm/hg19_rCRSchrm.fa";
my $interval = "/home/uec-00/shared/production/software/bissnp/wholegenome_interval_list.hg19.bed";
my $JAVA = "/home/uec-00/shared/production/software/java/default/bin/java -Xmx$ram" . "G";



my $cmd = "$JAVA -jar $bissnp -aecm -R $ref -T BisulfiteGenotyper -I $input -D $rod -vfn1 $input.cpg.raw.vcf -vfn2 $input.snp.raw.vcf -stand_call_conf 30 -stand_emit_conf 0 -L $interval -out_modes DEFAULT_FOR_TCGA -single_sample normal_test -nt $numcores -pem";

runcmd("$cmd");

sub runcmd
{
        my $cmd = shift @_;
        print STDERR "$cmd\n";
        system($cmd);
}
