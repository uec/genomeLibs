#!/usr/bin/perl
use File::Basename;
use lib dirname (__FILE__);
use EpigenomeUtils;

$outputPrefix = shift @ARGV;

$ENV{'BOWTIE_INDEXES'} = "$GENOMEROOT/genomes/bowtie/";
$ENV{'PATH'} .= "$SOFTWAREROOT/cufflinks/default:$SOFTWAREROOT/tophat/default:$SOFTWAREROOT/bowtie/default";

$execmd = "$SOFTWAREROOT/cufflinks/default/cuffdiff " . join(" ", @ARGV);
print "$execmd\n";
system($execmd);


@outputFiles = ("isoforms.fpkm_tracking", "genes.fpkm_tracking", "cds.fpkm_tracking","tss_groups.fpkm_tracking", "isoform_exp.diff","gene_exp.diff","tss_group_exp.diff","cds_exp.fpkm_tracking","splicing.diff","cds.diff","cds_exp.diff","promoters.diff");

for $f (@outputFiles)
{
	system("mv $f " . $outputPrefix. ".cuffdiff_$f");
}
