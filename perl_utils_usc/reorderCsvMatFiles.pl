#!/usr/bin/env perl

use strict;

$::DELIM = ",";
$::DEFAULT_VAL = "nan";

my $USAGE = "reorderCsvMatFiles.pl numHeaderCols file1.csv file2.csv ...";
die "$USAGE\n" unless (@ARGV>=2);
my ($numHeads, @fns) = @ARGV;

# Make a first pass to capture all ids
my $ids = {};
foreach my $fn (@fns)
{
    die "Can't read $fn\n" unless (open(F,$fn));
    my $thisIds = {};
    LINE: while (my $line = <F>)
    {
	next LINE if ($line =~ /^chrom/i); # header

	chomp $line;
	my @f = split(/$::DELIM/,$line);
	@f[0..($numHeads-1)] = map(&standardChrFld, @f[0..($numHeads-1)]);
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
    my $lineNum = 0;
    LINE: while (my $line = <F>)
    {
	next LINE if ($line =~ /^chrom/i); # header
	print STDERR "$fn\tline $lineNum\n" if (($lineNum % 10000)==0);
	chomp $line;
	my @f = split(/$::DELIM/,$line);
	@f[0..($numHeads-1)] = map(&standardChrFld, @f[0..($numHeads-1)]);
	my $id = fldsToId(\@f, $numHeads);
	#print STDERR "ID=$id\n";
	$lines->{$id} = join(",",@f); #\@f;  # Storing the array took lots of memory.
	
	print STDERR "Found two different column counts in $fn: " . $nFlds . ", " . scalar(@f) . "\n" if ( $nFlds && ($nFlds != scalar(@f)) );
	$nFlds = scalar(@f);
	$lineNum++;
    }
    close(F);

    # Now output in correct order
    my $outfn = $fn . ".ordered";
    die "Can't write to $outfn\n" unless (open(F,">$outfn"));
    foreach my $id (@orderedIds)
    {
	my $line = $lines->{$id};
	if (!$line)
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
    # Chrom fields have already been standardized
    # my $id = join($::DELIM, map (&standardChrFld,@$flds[0..($numHeads-1)]));
    my $id = join($::DELIM, @$flds[0..($numHeads-1)]);
    return $id;
}

sub standardChrFld
{
    my $fld = $_;
    my $outFld = $fld;
    if ($fld =~ /^chr(.*)$/)
    {
	$outFld = $1;
    }

    $outFld =~ s/M/23/gi;
    $outFld =~ s/X/24/gi;
    $outFld =~ s/Y/26/gi;
	
    return $outFld;
}
