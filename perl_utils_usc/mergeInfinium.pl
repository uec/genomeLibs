#!/usr/bin/env perl

use strict;

# Constants
my $WINLINE = chr(0xd);


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
	$line =~ s/$WINLINE//g;

	# Deal with excel quoted fields
	$line = escapeExcelSecs($line, $delim);

	@flds = split(/$delim/,$line);
#	@flds = map {$a=$_; $a=~s/^\"//g; $a=~s/\"$//g; $a} @flds;
	@flds = map {$a=$_; $a=$EMPTY_STR if ($a=~/^NA$/i); $a} @flds;
	my $masterLineNum = 0;
	if ($curLineNum == 1)
	{
	    # Header
	    $masterLineNum = 1;
	}
	else
	{
	    die "Found illegal line in $fn (must start with cg000000): $line\n" unless (($flds[0]=~/^(\d+[^,]*)/) || ($flds[0]=~/^(ch\.[^,]*)/) || ($flds[0]=~/^(rs\d[^,]*)/) || ($flds[0] =~ /^(cg\d+)/));
	    my $id = $1;
	    $lineNumsById->{$id} = $curLineNum if ($master);
	    
	    $masterLineNum = $lineNumsById->{$id};
	    if (!$masterLineNum)
	    {
#	    die "Can't find ID \"$id\" in master\n";
		next LINE;
	    }
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

    my $out = join(",",@{$lineCols})."\n";
    $out = unescapeExcelSecs($out,",");
    print STDOUT $out;
}



sub escapeExcelSecs($line)
{
    my ($line, $delim) = @_;

    my $rest = $line;
    my $out = "";
    my $nhits = 0;
    while ($rest =~ /^([^\"]*)\"([^\"]*)\"(.*)$/)
    {
	my ($a,$b,$c) = ($1,$2,$3);
	$b =~ s/${delim}/___/g;
	$out .= $a . $b;
	$rest = $c;
	$nhits++;
    }
    $out .= $rest;
#    print STDERR "GENERATED ESCAPED with ${nhits} matches: $out\n";

    return $out;
}

sub unescapeExcelSecs()
{
    my ($line, $delim) = @_;

#    print STDERR "UNESCAPING: $line\n";

    my @f = split($delim, $line);
    @f = map { my $out=$_; if ($out=~/___/) { $out =~ s/___/${delim}/g; $out = "\"${out}\"";  } $out; } @f;

    my $out = join($delim, @f);
#    print STDERR "\nUNESCAPED: $out\n";
    return $out;
}
