#!/usr/bin/perl

use strict;

my $USAGE = "fastaGffSubseqs.pl locs.gff chr1.fa chr2.fa ...";

die "$USAGE\n" unless (@ARGV>=2);
my ($gff_fns, @fa_fns) = @ARGV;


my $fa_map = readFastaSeqs(@fa_fns);

# Go through gff
die "Can't read $gff_fns\n" unless open(GFF,$gff_fns);
while (my $line = <GFF>)
{
    chomp $line;
    my @flds = split(/\t/,$line);
    my $chr = $flds[0];
    my ($s,$e) = @flds[3..4];
    my $strand = $flds[6];

#    print STDERR "Getting ${chr}:$s-$e ($strand)...";
    my $seq = getFastaSeq($fa_map, $chr, $s, $e, $strand);
#    print STDERR " $seq\n";

    my $ctcf_dir = ctcfDir($seq);
    $flds[6] = $ctcf_dir;

#    $flds[8] = $seq;
    print join("\t",@flds)."\n";
}
close (GFF);

sub ctcfDir
{
    my ($inseq) = @_;
    
    my @ev = (0,0);
    foreach my $ind (0,1)
    {
	my $seq = ($ind) ? $inseq : revcomp($inseq);

#	$ev[$ind]++ if (uc(substr($seq,3,1)) eq 'C');
	$ev[$ind]+=2 if (uc(substr($seq,4,2)) eq "CA");

	$ev[$ind]+=2 if (uc(substr($seq,8,2)) eq "AG");

#	$ev[$ind]++ if (uc(substr($seq,12,1)) eq 'G');

	$ev[$ind]+=2 if (uc(substr($seq,14,2)) eq "CA");
    }

    my $diff = $ev[1] - $ev[0];

    my $out = '.';
    if ( $diff >= 4 )
    {
	$out = "+";
    }
    elsif ( $diff <= -4 )
    {
	$out = "-";
    }

    my $outseq = ($out eq '-') ? revcomp($inseq) : $inseq;

    print STDERR join("\t",$inseq, $outseq,$ev[0],$ev[1], abs($diff), $out)."\n";

    return $out;
}

sub getFastaSeq
{
    my ($fa_map, $chr, $s, $e, $strand) = @_;

    my $chr_seq = $fa_map->{$chr};
    my $out = 0;
    if (!$chr_seq)
    {
	print STDERR "Can't find sequence \"$chr_seq\"\n";
    }
    else
    {
	if ($e>=length($chr_seq))
	{
	    print STDERR "Requesting sequence $s-$e from $chr, which only has " . length($chr_seq) ." nucs\n";
	}
	else
	{
	    $out = substr($chr_seq, $s-1,$e-$s+1);
	    $out = revcomp($out) if ($strand eq '-');
	}
    }

    return $out;
}

sub revcomp
{
    my ($seq) = @_;

    # Stupid case hack
    $seq = reverse(lc($seq));
    $seq =~ s/g/C/g;
    $seq =~ s/a/T/g;
    $seq =~ s/t/A/g;
    $seq =~ s/c/G/g;
    
    return $seq;
}

sub readFastaSeqs
{
    my (@fa_fns) = @_;

    my $map = {};
    foreach my $fa_fn (@fa_fns)
    {
	print STDERR "Reading $fa_fn\n";
	die "Can't read $fa_fn\n" unless (open(FA,$fa_fn));

	my $cur_seq_name = "";
	while (my $line = <FA>)
	{
	    chomp $line;
	    if ($line =~ /^>(.*)/)
	    {
		$cur_seq_name = $1;
		die "Saw $cur_seq_name multiple times\n" if ($map->{$cur_seq_name});
		$map->{$cur_seq_name} = "";
		print STDERR "\tFound seq $cur_seq_name\n";
	    }
	    else
	    {
		$map->{$cur_seq_name} .= $line;
	    }
	}

	close(FA);
    }
    return $map;
}

