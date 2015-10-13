package EpigenomeUtils;
use strict;
use warnings;

use Exporter;
our @ISA = 'Exporter';
our @EXPORT = qw($SOFTWAREROOT $SAMTOOLS $PICARD $JAVA $GATKSNP $PICARDTMP $MAQ $PUBLICATIONDATA $TMP_DIR $GENOMEROOT &runcmd);

our $SOFTWAREROOT = "/primary/vari/software";
our $SAMTOOLS = "$SOFTWAREROOT/samtools/samtools";
our $PICARD = "$SOFTWAREROOT/picard/default";
our $JAVA = "$SOFTWAREROOT/java/default/bin/java";
our $GATKSNP = "$SOFTWAREROOT/GATK2/default/GenomeAnalysisTK.jar";
our $MAQ = "$SOFTWAREROOT/maq-0.7.1/maq";
our $PUBLICATIONDATA = $SOFTWAREROOT;
our $GENOMEROOT = "/primary/vari/genomicdata";
our $TMP_DIR = "/secondary/tmp";

sub runcmd
{
        my $cmd = shift @_;
        print STDERR "$cmd\n";
        system($cmd);
}