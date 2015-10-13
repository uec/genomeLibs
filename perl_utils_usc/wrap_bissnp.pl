#!/usr/bin/perl
use strict;
use File::Basename;
use lib dirname (__FILE__);
use EpigenomeUtils;

my $numcores = `cat /proc/cpuinfo | grep processor -c`;
chomp $numcores;
#$numcores = $numcores / 2 - 1;

#my $ram = 4 * $numcores;
my $ram = 12;

my $input = $ARGV[0] || die "need input file";
my $bissnp = "$SOFTWAREROOT/bissnp/bissnp-default.jar";
my $vcf = "$SOFTWAREROOT/bissnp/dbsnp_135.hg19.sort.vcf";
my $ref = "GENOMEROOT/genomes/hg19_rCRSchrm/hg19_rCRSchrm.fa";
my $interval = "$SOFTWAREROOT/bissnp/wholegenome_interval_list.hg19.bed";
my $JAVA = "$SOFTWAREROOT/java/default/bin/java -Xmx$ram" . "G";



my $cmd = "$JAVA -jar $bissnp -R $ref -T BisulfiteGenotyper -I $input -D $vcf -vfn1 $input.cpg.raw.vcf -vfn2 $input.snp.raw.vcf -stand_call_conf 20 -stand_emit_conf 0 -L $interval -out_modes DEFAULT_FOR_TCGA -nt $numcores -rgv hg19 -mbq 0 -mmq 30";

runcmd("$cmd");

sub runcmd
{
        my $cmd = shift @_;
        print STDERR "$cmd\n";
        system($cmd);
}
