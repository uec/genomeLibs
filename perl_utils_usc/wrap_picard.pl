#!/usr/bin/perl
my $SAMTOOLS = "/home/uec-00/shared/production/software/samtools/samtools";
my $PICARD = "/home/uec-00/shared/production/software/picard/default/";
my $JAVA = "/home/uec-00/shared/production/software/java/default/bin/java";

my $cmd = join(" ", @ARGV);

runcmd("$JAVA -Xmx14g -jar $PICARD/$cmd");

sub runcmd
{
        my $cmd = shift @_;
        print STDERR "$cmd\n";
        system($cmd);
}