#!/usr/bin/env perl

use strict;


my $USAGE = "mergeById.pl file1.csv 1/2 file2.csv 5/6 ..."; # 1/2 means concatenate fields 1/2.  The last one does not have to be unique
my $DELIM = ",";


my $idFileCounts = {};
my $idToFieldHashes = []; # Each one is a hash for that particular file
my $nFiles = 0;
while (scalar(@ARGV)>0)
{
    my $fn = shift(@ARGV);
    my $idcol = shift(@ARGV);
    my $islast = (scalar(@ARGV)==0);

    print STDERR sprintf("file:%s, col:%d\n",$fn,$idcol);

    my $h = {};
    die "Can't read $fn\n" unless (open(F,$fn));
    LINE: while (my $line = <F>)
    {
	chomp $line;
	my @flds = split(/${DELIM}/,$line);
	my $id = concatFlds($idcol,\@flds);

	if (!$id)
	{
	    print STDERR "Why no id (col $idcol) for following line: $line\n";
	    next LINE;
	}

	if ($islast)
	{
	    # Output.  Only if we found the id in all other files.
	    if ($idFileCounts->{$id} < $nFiles)
	    {
		print STDERR "Why did we only find ID=$id a total of " . ($idFileCounts->{$id}) . " times? Skipping\n";
	    }
	    else
	    {
		my @newflds = (@flds);
		for (my $i=0; $i<$nFiles; $i++)
		{
		    my $otherh= @{$idToFieldHashes}[$i];
		    my $otherflds = $otherh->{$id};
		    @newflds = (@newflds,@$otherflds);
		}

		# Final output
		print join(",",@newflds)."\n";
	    }

	    if ($h->{$id})
	    {
	    }
	    else
	    {
		$idFileCounts->{$id}++;
	    }
	    
	}
	else
	{
	    # Store it
	    if ($h->{$id})
	    {
		print STDERR "Why duplicate ids ($id) in file $fn\n";
	    }
	    else
	    {
		$h->{$id} = \@flds;
		$idFileCounts->{$id}++;
	    }
	}

    }
    close(F);


    # Save and increment
    @$idToFieldHashes[$nFiles] = $h;
    $nFiles++;
}

sub concatFlds
{
    my ($concatColsSpec,$flds) = @_;

    my @cols = split(/\//,$concatColsSpec);
#    print STDERR "Cols = ".join(",",@cols)."\n";
    @cols = map {$_-1} @cols; # To make it an array index
    my $out = join("_",@{$flds}[@cols]);
    return $out;
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

