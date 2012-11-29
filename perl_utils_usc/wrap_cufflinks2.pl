#!/usr/bin/perl
use File::Basename;

#$ENV{'BOWTIE_INDEXES'} = "/home/uec-00/shared/production/genomes/bowtie/";
$ENV{'PATH'} .= "/home/uec-00/shared/production/software/cufflinks2/default:/home/uec-00/shared/production/software/tophat2/default:/home/uec-00/shared/production/software/bowtie2/default";
$file = $ARGV[$#ARGV];

$execmd = "/home/uec-00/shared/production/software/cufflinks2/default/cufflinks " . join(" ", @ARGV);
print "$execmd\n";
system($execmd);

$file = basename($file);
system("mv transcripts.gtf " . $file. ".cufflinks_transcripts.gtf");
system("mv transcripts.expr " . $file. ".cufflinks_transcripts.expr");
system("mv genes.expr " . $file. ".cufflinks_genes.expr");
