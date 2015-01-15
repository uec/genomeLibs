#!/usr/bin/perl
use File::Basename;
use lib dirname (__FILE__);
use EpigenomeUtils;

my $cmd = join(" ", @ARGV);

runcmd("$JAVA -Xmx12g -jar $PICARD/$cmd TMP_DIR=$PICARDTMP");

