#!/usr/bin/perl
use strict;

my $USAGE = "maqCompareOutput.pl full(0|1) minQ query.map.txt subject.map.txt";

die "$USAGE\n" if (@ARGV < 2);

my ($full_out, $minQ, $q_fn, $s_fn) = @ARGV;

my $all_hash = {};

my $q_hash = parseMaqMap($q_fn, $minQ);
print STDERR "query has " . (scalar(keys(%$q_hash))) . " entries\n";
foreach my $i (keys(%$q_hash)) { $all_hash->{$i} = 1 }


my $s_hash = parseMaqMap($s_fn, $minQ);
print STDERR "subject has " . (scalar(keys(%$s_hash))) . " entries\n";
foreach my $i (keys(%$s_hash)) { $all_hash->{$i} = 1 }

my $unmapped_q = 0;
my $unmapped_s = 0;
my $unmatched = 0;
foreach my $n (sort(keys(%$all_hash)))
{
    my $q = $q_hash->{$n};
    my $s = $s_hash->{$n};

    my @out = ($n);
    my $status = "neither";

    if ($q)
    {
	$status = "q_only";
	push(@out, @$q);
    }
    else
    {
	push (@out, map {'0'} (1..5));
	$unmapped_q++;
    }

    if ($s)
    {
	$status = "s_only";
	push(@out, @$s);
    }
    else
    {
	push (@out, map {'0'} (1..5));
	$unmapped_s++;
    }

    if ($q && $s)
    {

	my $diff = @{$q}[1] - @{$s}[1];
	push(@out, $diff);

	@{$q}[0] =~ /^chr(\d+)/;
	my $q_chr = $1;
	@{$s}[0] =~ /^chr(\d+)/;
	my $s_chr = $1;
	my $match = ($q_chr == $s_chr);

	$unmatched++ unless ($match);

	$status = ($match) ? "match" : "mismatch";
    }

    unshift(@out, $status);
    print join("\t",@out)."\n";
}

my $t =  scalar(keys(%$q_hash));
print "Unmapped_Q = $unmapped_q/$t (" . ($unmapped_q/$t) . ")\n";
print "Unmapped_S = $unmapped_s/$t (" . ($unmapped_s/$t) . ")\n";
print "Unmatched chroms = $unmatched\n";

# -- - - Funcs


# hash indexed by sequence name.  Each contains
# [$chr, $pos, $strand, $q, $c_or_t]
sub parseMaqMap
{
    my ($fn, $minQ) = @_;
    my $out = {};

    die "Can't read $fn\n" unless (open(F,$fn));

    my $line_n = 1;
    my $total_contam = 0;
  LINE: while (my $line = <F>)
  {
      chomp $line;
      my @flds = split(/\t/,$line);

      my $seq_name = @flds[0];
      my $locus = @flds[1];
      my $pos_raw = @flds[2];
      my $strand = @flds[3];
      my $seq = @flds[14];
      my $q = @flds[7];
      my $c_or_t = lc(substr($seq,0,1));

      next LINE if ($q < $minQ);

#    print "$locus\t$c_or_t\n";

      $locus =~ /(chr[0-9XYM]+.*)/;
      my $chr = $1;
      if (!$chr)
      {
	print STDERR "Bad chr: $locus\n";
	  $total_contam++;
	  next LINE;
      }
      my $key = $seq_name;

      my $pos = $pos_raw;
      if ($locus =~ /start\:(\d+)_/)
      {
	  $pos = $1;
      }
    

      print STDERR "Saw $key more than once\n" if ($out->{$key});
      $out->{$key} = [$chr, $pos, $strand, $q, $c_or_t];
  }

    close(F);

    return $out;
}

