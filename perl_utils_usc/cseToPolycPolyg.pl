#!/usr/bin/perl

use strict;
use DBI qw(:sql_types);

# Takes $numWinds windows on either side of the point.   Tabulates number of polyCs and polyGs in
# each window.
my $USAGE = "cseToPolycPolyg.pl 10 [windSize] 50 [numWinds] (cse file on STDIN)";
my $DBNAME = "genomes";
my $TABLE = "hg18genomeSeqs";

die "$USAGE\n" unless (@ARGV>=1);

my ($windSize, $numWinds) = @ARGV;

# DB connect
$::dbh = DBI->connect( "dbi:mysql:${DBNAME}",
'benb',
'',
{ AutoCommit => 0 }
    ) || die "Database connection not made: $DBI::errstr";


eval
{
    my $sql = "select avg(residue=?), count(*)  FROM ${TABLE} WHERE (chrom=?) AND (chromPos>=?) AND (chromPos<=?);";


    $::selectSth = $::dbh->prepare($sql);

    my $on_line = 1;
    my $balanceOfPower;
    while (my $line = <STDIN>)
    {
	print STDERR "On line $on_line (balance $balanceOfPower)\n" if (($on_line++ % 500)==0);
	chomp $line;
	my @flds = split(/[\t,]/,$line);
	if ((@flds < 2) || (@flds > 3))
	{
	    print STDERR "Line does not match c\ts\te format (C [\t,] S [\t,] [E])\n";
	    next;
	}
	
	my $chr = $flds[0];
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
	    my $count = getPolyCount($chr, $windStart, $windEnd);
	    push(@outFlds, $count);

	    # Figure out which half has more.  Right hand side is positive
	    my $balanceAdjust = $count * (($windStart<$center) ? -1 : 1);
	    $balanceOfPower += $balanceAdjust;
#	    print STDERR "$chr\t$windStart\t$windEnd\t$count\t$balanceOfPower\n";

	    # Increment
	    $windStart += $windSize;

	    # Skip the central CpG
	    # $windStart++ if ($windStart == $center);
	}
    
	# Reverse depending on bigger side
#	print STDERR "Balance:\t$balanceOfPower\n";
	if ($balanceOfPower < 0)
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



sub getPolyCount
{
    my ($chr, $s, $e) = @_;

    my $polycCount = getPolyCountHelper($chr, $s, $e, 'c');
    my $polygCount = getPolyCountHelper($chr, $s, $e, 'g');

    # Return the max
    return ($polycCount>$polygCount) ? $polycCount : $polygCount;
}

sub getPolyCountHelper
{
    my ($chr, $s, $e, $residue) = @_;

    $::selectSth->bind_param(1,$residue);
    $::selectSth->bind_param(2,$chr);
    $::selectSth->bind_param(3,$s);
    $::selectSth->bind_param(4,$e);

    $::selectSth->execute();

    my ($residueFrac, $numRows);
    $::selectSth->bind_columns(undef, \$residueFrac, \$numRows);
    $::selectSth->fetch();


    my $out = "NaN";
    if ($numRows>0)
    {
	$out = int(($residueFrac**2)*100); # Percentage
    }

#    print STDERR "$residue\t$chr\t$s\t$e\tResiudeFrac=${residueFrac}\tNumRows=$numRows\t$out\n";
    return $out;
}

