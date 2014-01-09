#!/usr/bin/perl

$ENV{PYTHONPATH} = "/home/uec-00/ramjan/lib/python2.7/site-packages";
$MACS = "/auto/uec-00/shared/production/software/MACS/1.4.2/bin/macs14";
$cmd = $MACS . " " . join(" ",@ARGV);

runcmd($cmd);

sub runcmd
{
        my $cmd = shift @_;
        print STDERR "$cmd\n";
        system($cmd);
}
