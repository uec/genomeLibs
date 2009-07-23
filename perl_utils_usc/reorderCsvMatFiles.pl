#!/usr/bin/env perl

use strict;

$::DELIM = ",";
$::DEFAULT_VAL = 0.0;

my $USAGE = "reorderCsvMatFiles.pl numHeaderCols file1.csv file2.csv ...";
die "$USAGE\n" unless (@ARGV>=2);
my ($numHeads, @fns) = @ARGV;

# Make a first pass to capture all ids
my $ids = {};
foreach my $fn (@fns)
{
    die "Can't read $fn\n" unless (open(F,$fn));
    my $thisIds = {};
    while (my $line = <F>)
    {
	chomp $line;
	my @f = split(/$::DELIM/,$line);
	my $id = fldsToId(\@f, $numHeads);
	# print STDERR "\tID=$id\n";
	$ids->{$id}++;
	$thisIds->{$id}++;

	print STDERR "Found ID multiple times in $fn: " . $id . "\n" if ($thisIds->{$id} == 2);
    }
    close(F);
    
    print STDERR "Read $fn,\t" .scalar(keys(%$thisIds)) . " in this file,\t" . scalar(keys(%$ids)) . " ids total\n";
}
my @orderedIds = sort(keys(%$ids));

# Second pass to output new files with correct ordering
foreach my $fn (@fns)
{
    # First index all lines
    die "Can't read $fn\n" unless (open(F,$fn));
    my $lines = {};
    my $nFlds = (); # Keep track so we can make empty lines
    while (my $line = <F>)
    {
	chomp $line;
	my @f = split(/$::DELIM/,$line);
	my $id = fldsToId(\@f, $numHeads);
	$lines->{$id} = \@f;
	
	print STDERR "Found two different column counts in $fn: " . $nFlds . ", " . scalar(@f) . "\n" if ( $nFlds && ($nFlds != scalar(@f)) );
	$nFlds = scalar(@f);
    }
    close(F);

    # Now output in correct order
    my $outfn = $fn . ".ordered";
    die "Can't write to $outfn\n" unless (open(F,">$outfn"));
    foreach my $id (@orderedIds)
    {
	my $f = $lines->{$id};
	my $line;
	if ($f)
	{
	    $line = join($::DELIM,@$f);
	}
	else
	{
	    # Make an empty line
	    $line = join($::DELIM, $id, map {$::DEFAULT_VAL} (1..($nFlds-$numHeads)));
	}

	print F $line . "\n";
    }
    close(F);
}




sub fldsToId
{
    my ($flds, $numHeads) = @_;
    my $id = join($::DELIM, @$flds[0..($numHeads-1)]);
    return $id;
}
