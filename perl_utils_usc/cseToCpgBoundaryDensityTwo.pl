#!/usr/bin/perl

use strict;
use DBI qw(:sql_types);

# Takes $numWinds windows on either side of the CpG.   Tabulates number of Cpgs in windows.
my $USAGE = "cseToCpgBoundaryDensity.pl 10 [windSize] 50 [numWinds] 1 [auto flip] < (cse file on STDIN)\nauto-flip, if set, will flip so that the higher density half is on the right";
my $DBNAME = "genomes";
my $TABLE = "hg18genomeSeqs";

die "$USAGE\n" unless (@ARGV>=1);

my ($windSize, $numWinds, $autoFlip) = @ARGV;

# DB connect
$::dbh = DBI->connect( "dbi:mysql:${DBNAME}",
'benb',
'',
{ AutoCommit => 0 }
    ) || die "Database connection not made: $DBI::errstr";


eval
{
    my $sql = "SELECT avg(isCpG), count(isCpG) FROM ${TABLE} WHERE (chrom=?) AND (chromPos>=?) AND (chromPos<=?);";
#    my $sql = "SELECT GROUP_CONCAT(isCpG) FROM ${TABLE} WHERE (chrom=?) AND (chromPos>=?) AND (chromPos<=?);";

    $::selectSth = $::dbh->prepare($sql);

    my $on_line = 1;
    my $balanceOfPower;
    while (my $line = <STDIN>)
    {
	print STDERR "On line $on_line (balance $balanceOfPower)\n" if (($on_line++ % 100)==0);
	chomp $line;
	my @flds = split(/[\t,]/,$line);
	if ((@flds < 2) || (@flds > 3))
	{
	    print STDERR "Line does not match c\ts\te format (C [\t,] S [\t,] [E])\n";
	    next;
	}
	print STDERR join(",",@flds[0..3])."\n";
	
	my $chr = $flds[0];
	if ($chr == 23)
	{
		$chr = 'M';	
	}
	elsif ($chr == 24)
	{
		$chr = 'X';	
	}
	elsif ($chr == 26)
	{
		$chr = 'Y';	
	}
	$chr = "chr$chr" unless ($chr =~ /^chr/);
		
	my $center = $flds[1];
	my @outFlds = ();

	my $first = $center - ($windSize*$numWinds);
	my $last = $center + ($windSize*$numWinds);

	my $windStart = $first;
	$balanceOfPower = 0;
	while ($windStart < $last)
	{
	    # Get CpG count
	    my $windEnd = $windStart + $windSize - 1;
	    my $count = getCpgCount($chr, $windStart, $windEnd);
	    push(@outFlds, $count);

	    # Figure out which half has more.  Right hand side is positive
	    my $balanceAdjust = $count * (($windStart<$center) ? -1 : 1);
	    $balanceOfPower += $balanceAdjust;
	    #print STDERR "$chr\t$windStart\t$windEnd\t$count\t$balanceOfPower\n";

	    # Increment
	    $windStart += $windSize;
	    # Skip the central CpG
	    $windStart++ if ($windStart == $center);
	}
    
	# Reverse depending on bigger side
#	print STDERR "Balance:\t$balanceOfPower\n";
	if ($autoFlip && ($balanceOfPower < 0))
	{
	    @outFlds = reverse(@outFlds);
	}

	print join(",",@outFlds)."\n";

    }

};
if( $@ ) {
    print STDERR "Database error: $DBI::errstr\n";
    $::dbh->rollback(); #just die if rollback is failing
}


# Disconnect from DB
$::selectSth->finish();
$::dbh->disconnect();



sub getCpgCount
{
    my ($chr, $s, $e) = @_;

    # The CpG is always on the G, so shift forward by one
    $s++;
    $e++;

	# For some reason, DBI chokes with negative start/ends
	my $numRows = 0;
	my $cpgCount;
	if (($s>0) && ($e>0))
	{ 

    $::selectSth->bind_param(1,$chr);
    $::selectSth->bind_param(2,$s);
    $::selectSth->bind_param(3,$e);

   $::selectSth->execute();

    $::selectSth->bind_columns(undef, \$cpgCount, \$numRows);
    $::selectSth->fetch();
#   print STDERR "$chr\t$s\t$e\tCpgCount=${cpgCount}\tNumRows=$numRows\n";
	}
	
    return ($numRows>0) ? int($cpgCount*100) : "NaN";  # Percentage
}

