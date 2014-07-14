#!/usr/bin/perl
use File::Basename;

$ENV{'BOWTIE_INDEXES'} = "/home/uec-00/shared/production/genomes/bowtie/";
$ENV{'PATH'} .= "/home/uec-00/shared/production/software/cufflinks/default:/home/uec-00/shared/production/software/tophat/default:/home/uec-00/shared/production/software/bowtie/default";

$execmd = "/home/uec-00/shared/production/software/cufflinks/default/cuffcompare " . join(" ", @ARGV);
print "$execmd\n";
system($execmd);
