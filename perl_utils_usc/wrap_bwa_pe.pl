#!/usr/bin/perl -w
# Mon Jul 12 14:03:06 PDT 2010
use strict;

# - assumes paired end data, see wrap_bwa_se.pl for single end alignment
# - assumes SA file will be named as read1.sai from read1.fq
my $usage = "wrap_bwa_pe.pl refFa.fa read1.fq read2.fq [outfile name]";
die "$usage\n" unless (@ARGV >= 3);

my $bwa = "/auto/uec-00/shared/production/software/bwa/default/bwa";
#my $bwa = "/auto/uec-01/shared/knowles/software_src/bwa-0.5.8a/bwa";

my $refFa = $ARGV[0];
my $read1 = $ARGV[1];
my $read2 = $ARGV[2];

my $read1SA = $read1;
my $read2SA = $read2;

$read1SA =~ s/\.fq$//i;
$read2SA =~ s/\.fq$//i;

my $outfileSAM = $ARGV[3] || $read1SA . "." . $read2SA . ".mapped.sam";

if ($outfileSAM !~ /.sam$/)
{
	$outfileSAM .= ".sam";
}

$read1SA .= ".sai";
$read2SA .= ".sai";

# check existence of ref genome
if ( ! -e $refFa )
{
	die "reference does not exist.\n";
}

# check existence of read sequence files
die "need read sequence files\n" unless ( -e $read1 && -e $read2 );

# when calling, redirect STDERR?
my $cmd = join(" ", $bwa, "aln", "-I -t 8", $refFa, $read1, "> $read1SA");
# system call
#print "$cmd\n";
system($cmd);

$cmd = join(" ", $bwa, "aln", "-I -t 8", $refFa, $read2, "> $read2SA");
# system call
#print "$cmd\n";
system($cmd);

# check existence of SA files
die "need read.sai files\n" unless ( -e $read1SA && -e $read2SA );

my $alignCMD = join(" ", $bwa, "sampe", $refFa, $read1SA, $read2SA, $read1, $read2, "> $outfileSAM");
# system call
#print "$alignCMD\n";
system($alignCMD);

