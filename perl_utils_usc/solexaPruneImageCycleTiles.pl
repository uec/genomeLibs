#!/usr/bin/perl

use strict;
use File::Basename qw/basename/;

my $USAGE = "solexaPruneImageCycleTiles [-delete | -test] RunDir/Data/Images cycles-mod tiles-mod";

my $confirm = $ARGV[0] || die "must run with either -delete or -test flag\n$USAGE";
$confirm =~ /\-test/ || $confirm =~ /\-delete/ ||die "must run with either -delete or -test flag\n$USAGE";

my $runDir = $ARGV[1] || die "specify dir to run on: ex /srv/data/slxa/incoming/myrun123AAXX/Images\n$USAGE";
-d $runDir || die "DIR NOT FOUND. specify dir to run on: ex /srv/data/slxa/incoming/myrun/Images\n$USAGE";

my $cycleMod = $ARGV[2] || die "specify cycles modulus to keep\n$USAGE";
$cycleMod =~ /\d+/ || die "specify cycles modulus to keep\n$USAGE";

my $tileMod = $ARGV[3] || die "specify tile modulus to keep\n$USAGE";
$tileMod =~ /\d+/ || die "specify tile modulus to keep\n$USAGE";

my @files = glob("$runDir/L00?/C*/*.tif");

foreach my $tif (@files)
{
	$tif =~ /L00(\d)\/C(\d+)\.\d\/s_\d_(\d+)_\w\.tif/;
	my $laneNumber = $1;
	my $cycleNumber = $2;
	my $tileNumber = $3;
	if(($cycleNumber-1) % $cycleMod != 0 || $tileNumber % $tileMod != 0)
	{
		if($confirm =~ /\-delete/)
		{
			print "deleting $tif\n";
			unlink($tif);
		}
		else 
		{
#			print "would have deleted $tif\n";
		}
	}
	else
	{
		print "keeping  $tif\n";
	}
}

exit;
