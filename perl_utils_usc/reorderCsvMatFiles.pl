#!/usr/bin/env perl

use strict;

$::DELIM = ",";
$::DEFAULT_VAL = "null";

my $USAGE = "reorderCsvMatFiles.pl numHeaderCols file1.csv file2.csv ...";
die "$USAGE\n" unless (@ARGV>=2);
my ($numHeads, @fns) = @ARGV;

# Make a first pass to capture all ids
my $ids = {};
foreach my $fn (@fns)
{
    die "Can't read $fn\n" unless (open(F,$fn));
    my $thisIds = {};
    my $onLine = 0;
    LINE: while (my $line = <F>)
    {
	#next LINE if (($line =~ /^chrom/i) || ($line=~/ample/)); # header
	$onLine++;
	next LINE if ($onLine==1);

	chomp $line;
	my @f = split(/$::DELIM/,$line);
	@f[0..($numHeads-1)] = map(&standardChrFld, @f[0..($numHeads-1)]);
	my $id = fldsToId(\@f, $numHeads);
	next LINE if ($id =~ /^\s*$/);
	
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
    my $headerLine = "";
    LINE: while (my $line = <F>)
    {
	$lineNum++;
	
	print STDERR "$fn\tline $lineNum\n" if (($lineNum % 10000)==0);
	chomp $line;
	
	# Get rid of windows line feeds
	my @lineParts = split(chr(0xd),$line);
	$line = $lineParts[0];
	my @f = split(/$::DELIM/,$line);
	my $newLine = join(",",@f);
	
	if ($lineNum == 1)
	{
		$headerLine = $newLine;
		next LINE;
	}
	
	
	
	@f[0..($numHeads-1)] = map(&standardChrFld, @f[0..($numHeads-1)]);
	my $id = fldsToId(\@f, $numHeads);
	#print STDERR "ID=$id\n";
	
	# Did we see this ID before?
	my $oldLine = $lines->{$id};
	my $lineToUse = ($oldLine) ? betterLine($oldLine, $newLine) : $newLine;
	
	$lines->{$id} = $lineToUse; #\@f;  # Storing the array took lots of memory.
	
	print STDERR "Found two different column counts in $fn: " . $nFlds . ", " . scalar(@f) . "\n" if ( $nFlds && ($nFlds != scalar(@f)) );
	$nFlds = scalar(@f);
    }
    close(F);

    # Now output in correct order
    my $outfn = $fn . ".ordered";
    die "Can't write to $outfn\n" unless (open(F,">$outfn"));
    print F $headerLine ."\n";
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
    
    my @stdFlds = map { standardHeadFld($_); } @$flds;
    
    my $id = join($::DELIM, @stdFlds[0..($numHeads-1)]);
    return $id;
}

sub standardHeadFld
{
	my ($inFld) = @_;
	
	my $outFld = $inFld;
	
	# Special case
	if ($outFld =~ /TCGA/)
	{
		$outFld =~ s/_/-/g;
		$outFld = uc($outFld);
		if ($outFld !~ /(TCGA-\w\w+-\w\w\w\w+)/)
		{
			print STDERR "Illegal TCGA fomat: $outFld\n";
		}
		else
		{
			$outFld = $1;
		}
	}

	return $outFld;
}

sub betterLine
{
	my ($oldLine, $newLine) = @_;
	my $lineToUse = $oldLine; # By default, use the first one.
	
	my @oldFlds = split(/,/,$oldLine);
	my @newFlds = split(/,/,$newLine);
	
	# Special case, TCGA slide file
	if (($oldFlds[14] eq "TOP") || ($oldFlds[14] eq "BOTTOM"))
	{
		print STDERR sprintf("%s: Comparing %d to %d\n",$oldFlds[0],$oldFlds[13],$newFlds[13]);
		
		# Use the lower of the two
		$lineToUse = ($oldFlds[13] <= $newFlds[13]) ? $oldLine : $newLine;
	}
	
	return $lineToUse;
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
