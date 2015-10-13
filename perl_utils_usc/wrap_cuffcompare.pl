#!/usr/bin/perl
use File::Basename;
use lib dirname (__FILE__);
use EpigenomeUtils;

$ENV{'BOWTIE_INDEXES'} = "$GENOMEROOT/genomes/bowtie/";
$ENV{'PATH'} .= "$SOFTWAREROOT/cufflinks/default:$SOFTWAREROOT/tophat/default:$SOFTWAREROOT/bowtie/default";

$execmd = "$SOFTWAREROOT/cufflinks/default/cuffcompare " . join(" ", @ARGV);
print "$execmd\n";
system($execmd);
