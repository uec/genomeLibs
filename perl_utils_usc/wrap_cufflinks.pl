#!/usr/bin/perl
use File::Basename;
use lib dirname (__FILE__);
use EpigenomeUtils;

$ENV{'BOWTIE_INDEXES'} = "$GENOMEROOT/genomes/bowtie/";
$ENV{'PATH'} .= "$SOFTWAREROOT/cufflinks/default:$SOFTWAREROOT/tophat/default:$SOFTWAREROOT/bowtie/default";
$file = $ARGV[$#ARGV];

$execmd = "$SOFTWAREROOT/cufflinks/default/cufflinks " . join(" ", @ARGV);
print "$execmd\n";
system($execmd);

$file = basename($file);
system("mv transcripts.gtf " . $file. ".cufflinks_transcripts.gtf");
system("mv transcripts.expr " . $file. ".cufflinks_transcripts.expr");
system("mv genes.expr " . $file. ".cufflinks_genes.expr");
