#!/usr/bin/perl
use strict;
use File::Basename;
use lib dirname (__FILE__);
use EpigenomeUtils;

$ENV{PYTHONPATH} = "/home/uec-00/ramjan/lib/python2.7/site-packages";
my $MACS = "$SOFTWAREROOT/MACS/1.4.2/bin/macs14";
my $cmd = $MACS . " " . join(" ",@ARGV);

runcmd($cmd);
