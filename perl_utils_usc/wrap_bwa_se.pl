#!/usr/bin/perl -w
use strict;

# - assumes single end data, see wrap_bwa_pe.pl for paired end alignment
# - assumes SA file will be named as read1.sai from read1.fq
my $USAGE = "wrap_bwa_se.pl refFa.fa read1.fq [outfile name]";
die "$USAGE\n" unless (@ARGV >= 2);

my $bwa = "/auto/uec-00/shared/production/software/bwa/default/bwa";

my $refFa = $ARGV[0];
my $read1 = $ARGV[1];

my $read1SA = $read1;

$read1SA =~ s/\.fq$//i;

my $outfileSAM = $ARGV[2] || $read1SA . ".mapped.sam";

if ($outfileSAM !~ /.sam$/)
{
	$outfileSAM .= ".sam";
}

$read1SA .= ".sai";

# check existence of ref genome
if ( ! -e $refFa )
{
	die "reference does not exist.\n";
}

# check existence of read sequence files
die "need read sequence files\n" unless ( -e $read1 );

# when calling, redirect STDERR?
my $cmd = join(" ", $bwa, "aln", "-I -t 8", $refFa, $read1, "> $read1SA");
print "$cmd\n";
# system call
system($cmd);

my $alignCMD = join(" ", $bwa, "samse", $refFa, $read1SA, $read1, "> $outfileSAM");
print "$alignCMD\n";
# system call
system($alignCMD);


