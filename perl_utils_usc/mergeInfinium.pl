#!/usr/bin/env perl

use strict;


my @fns = (@ARGV);

my $EMPTY_STR = "nan";

my $lineNumsById = 0;
my $outLines = [];
my $curStartCol = 0;
foreach my $fn (@fns)
{
    my $delim = ($fn =~ /.txt$/) ? "\t" : ",";

    my $master = (!$lineNumsById);
    $lineNumsById = {} if ($master);
    
    die "Can't read $fn\n" unless (open(F,$fn));
    my $curLineNum = 1;
    my @flds = ();
    my $cgsSeen = {};
    LINE: while (my $line = <F>)
    {
	chomp $line;
	@flds = split(/$delim/,$line);
	my $masterLineNum = 0;
	if ($curLineNum == 1)
	{
	    # Header
	    $masterLineNum = 1;
	}
	else
	{
	    die "Found illegal line in $fn (must start with cg000000): $line\n" unless ($line =~ /^(cg\d+)/);
	    my $id = $1;
	    $lineNumsById->{$id} = $curLineNum if ($master);
	    
	    $masterLineNum = $lineNumsById->{$id};
	    die "Can't find ID \"$id\" in master\n" unless ($masterLineNum);
	    next LINE if ($cgsSeen->{$id});  # Houtan's spreadsheet had some duplicate lines.
	    $cgsSeen->{$id} = 1;
	}

	# Fill in any empty fields
	@{$outLines}[$masterLineNum-1] = [] if ($master);
	my $colsSoFar = scalar(@{$outLines}[$masterLineNum-1]);
	for (my $i = $colsSoFar; $i < $curStartCol; $i++)
	{
	    @{@${outLines}[$masterLineNum-1]}[$i] = $EMPTY_STR;
	}

	# Add the new fields
	push(@{@${outLines}[$masterLineNum-1]}, @flds);

	$curLineNum++;
    }

    $curStartCol += scalar(@flds);

    close(F);
}



foreach my $lineCols (@${outLines})
{
# Go through and complete the last file's columns
    my $colsSoFar = scalar(@{$lineCols});
    for (my $i = $colsSoFar; $i < $curStartCol; $i++)
    {
	@{$lineCols}[$i] = $EMPTY_STR;
    }

    print STDOUT join(",",@{$lineCols})."\n";
}
