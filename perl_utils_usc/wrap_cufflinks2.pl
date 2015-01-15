#!/usr/bin/perl
use File::Basename;

#$ENV{'BOWTIE_INDEXES'} = "/home/uec-00/shared/production/genomes/bowtie/";
$ENV{'PATH'} .= "$SOFTWAREROOT/cufflinks2/default:$SOFTWAREROOT/tophat2/default:$SOFTWAREROOT/bowtie2/default";
$file = $ARGV[$#ARGV];

$execmd = "$SOFTWAREROOT/cufflinks2/default/cufflinks " . join(" ", @ARGV);
print "$execmd\n";
system($execmd);

$file = basename($file);
system("mv transcripts.gtf " . $file. ".cufflinks_transcripts.gtf");
system("mv isoforms.fpkm_tracking " . $file. ".cufflinks_isoforms.expr");
system("mv genes.fpkm_tracking " . $file. ".cufflinks_genes.expr");
