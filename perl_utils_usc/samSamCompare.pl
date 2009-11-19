#!/usr/bin/env perl

use strict;

my $USAGE = "samSamCompare.pl MINQ file1.sam file2.sam";

die "$USAGE\n" unless (@ARGV == 3);
my ($minq, $fn1, $fn2) = @ARGV;

print STDERR "Reading $fn1...\n";
my $sams1 = mapSamReads($fn1, $minq);
print STDERR "Reading $fn2...\n";
my $sams2 = mapSamReads($fn2, $minq);



my $readsSeen = {};
print STDERR "Processing $fn1 alignments...\n";
foreach my $key (sort keys(%$sams1))
{
    my $sam1 = $sams1->{$key};
    my $sam2 = $sams2->{$key};

    my $code = ".";
    if ($sam2)
    {
	$code = samsEqual($sam1,$sam2) ? "match" : "mismatch";
    }
    else
    {
	$code = "missing2";
	$sam2 = [0, 0, 0];
    }

    print join("\t",$key,$code);
    print "\t".join("\t",@$sam1,@$sam2) unless ($code eq 'match');
    print "\n";

    $readsSeen->{$key}++;
}


print STDERR "Processing $fn2 alignments...\n";
KEY2: foreach my $key (sort keys(%$sams2))
{
    my $seen = $readsSeen->{$key}++;
    next KEY2 if ($seen);

    my $code = "missing1";
    my $sam1 = [0, 0, 0];
    my $sam2 = $sams2->{$key};
    print join("\t",$key,$code);
    print "\t".join("\t",@$sam1,@$sam2);
    print "\n";
}



# map values are lists [chrom, rev-strand(0|1), startCoord]
sub mapSamReads
{
    my ($fn, $minq) = @_;
    my $out = {};

    die "Can't read $fn\n" unless (open(F,$fn));
    LINE: while (my $line = <F>)
    {
	chomp $line;
	my @fld = split("\t",$line);

	next LINE if ($fld[4] < $minq);
	
	my $key = $fld[0];
#	$key =~ s/\/[12]$//g;


	my $sam = [ $fld[2],
		    (($fld[1] & 16)>0) ? 1:0,
		    $fld[3] ];
	$out->{$key} = $sam;
    }
        
    close(F);
    return $out;
}

sub samsEqual
{
    my ($sam1, $sam2) = @_;

    return 0 unless (@$sam1[0] eq @$sam2[0]);
    return 0 unless (@$sam1[1] == @$sam2[1]);
    return 0 unless (@$sam1[2] == @$sam2[2]);
    return 1;
}

