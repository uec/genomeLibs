#!/usr/bin/env perl

use strict;
use Getopt::Long;

my $USAGE = "bedToUcscHyperlinks.pl [--prefix http://epiweb.usc.edu/] [--genome hg18] < a.bed > a.txt";

my $prefix = "http://epiweb.usc.edu/";
my $genome = "hg18";
GetOptions ('prefix=s' => \$prefix, 'genome=s' => \$genome) || print STDERR "$USAGE\n";
print STDERR "prefix=$prefix\tgenome=$genome\n";

# Headers
my @headers = ("UCSC link","chr","start","end","name","score","strand");
print join("\t",@headers)."\n";


my $itemCount = 0;
LINE: while (my $line = <>)
{
    chomp $line;
    my @flds = split(/\t/,$line);
    my ($chr, $s, $e, $name, $score, $strand,@rest) = @flds;

    if (!isInt($s) || !isInt($e))
    {
	print STDERR "Skipping non-BED line: $line\n";
	next LINE;
    }

    $itemCount++;

    # Sometimes only use 4 columns and put the score as the last column.
    # Sloppy BED, but we try to detect if
    if (isNumeric($name) && !$score)
    {
	$score = $name;
	$name = '';
    }

    if (!$name)
    {
	$name = "el${itemCount}";
    }

    my $url = makeUrl($prefix, $genome, $chr, $s, $e);
    my $hl = makeExcelHyperlink($url, $name);

    print join("\t",($hl, $chr,$s,$e,$name,$score,$strand,@rest))."\n"
}

sub makeExcelHyperlink
{
    my ($url,$name) = @_;

    my $hl = "=HYPERLINK(\"$url\",\"$name UCSC link\")";
}

sub makeUrl
{
    my ($prefix, $genome, $chr, $start, $end) = @_;

    my $url = "${prefix}cgi-bin/hgTracks?position=${chr}:${start}-${end}&db=${genome}";  
}

sub isInt
{
    my ($var) = @_;
    
    return ($var =~ /^\d+$/);
}

sub isNumeric
{
    my ($var) = @_;
    
    return ($var =~ /^[0-9\.\-E]+$/);
}
