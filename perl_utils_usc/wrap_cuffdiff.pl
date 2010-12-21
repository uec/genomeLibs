#!/usr/bin/perl
use File::Basename;

$outputPrefix = shift @ARGV;

$ENV{'BOWTIE_INDEXES'} = "/home/uec-00/shared/production/genomes/bowtie/";
$ENV{'PATH'} .= "/home/uec-00/shared/production/software/cufflinks/default:/home/uec-00/shared/production/software/tophat/default:/home/uec-00/shared/production/software/bowtie/default";

$execmd = "/home/uec-00/shared/production/software/cufflinks/default/cuffdiff " . join(" ", @ARGV);
print "$execmd\n";
system($execmd);


@outputFiles = ("isoforms.fpkm_tracking", "genes.fpkm_tracking", "cds.fpkm_tracking","tss_groups.fpkm_tracking", "isoform_exp.diff","gene_exp.diff","tss_group_exp.diff","cds_exp.fpkm_tracking","splicing.diff","cds.diff","cds_exp.diff","promoters.diff");

for $f (@outputFiles)
{
	system("mv $f " . $outputPrefix. ".cuffdiff_$f");
}
