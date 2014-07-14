#!/usr/bin/perl
my $SAMTOOLS = "/home/uec-00/shared/production/software/samtools/samtools";
my $PICARD = "/home/uec-00/shared/production/software/picard/default";
my $JAVA = "/home/uec-00/shared/production/software/java/default/bin/java";

my $cmd = join(" ", @ARGV);

runcmd("$JAVA -Xmx12g -jar $PICARD/$cmd TMP_DIR=/export/uec-gs1/laird/shared/tmp");

sub runcmd
{
        my $cmd = shift @_;
        print STDERR "$cmd\n";
        system($cmd);
}
