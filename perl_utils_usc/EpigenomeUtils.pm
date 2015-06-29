package EpigenomeUtils;
use strict;
use warnings;

use Exporter;
our @ISA = 'Exporter';
our @EXPORT = qw($SOFTWAREROOT $SAMTOOLS $PICARD $JAVA $GATKSNP $PICARDTMP $MAQ $PUBLICATIONDATA $TMP_DIR &runcmd);

our $SOFTWAREROOT = "/home/uec-00/shared/production/software";
our $SAMTOOLS = "$SOFTWAREROOT/samtools/samtools";
our $PICARD = "$SOFTWAREROOT/picard/default";
our $JAVA = "$SOFTWAREROOT/java/default/bin/java";
our $GATKSNP = "$SOFTWAREROOT/GATK2/default/GenomeAnalysisTK.jar";
our $MAQ = "$SOFTWAREROOT/maq-0.7.1/maq";
our $PUBLICATIONDATA = "/home/uec-00/shared/publicationData";
our $TMP_DIR = "/secondary/tmp";

sub runcmd
{
        my $cmd = shift @_;
        print STDERR "$cmd\n";
        system($cmd);
}