#!/usr/bin/perl

use strict;
use File::Basename qw/basename/;
use threads;
use Thread::Semaphore;

my $semaphore = Thread::Semaphore->new(8);
my $workDir = "slxaMagickTmp";

my $USAGE = "solexaMagick RunDir/Data/Images ";

my $runDir = $ARGV[0] || die "specify dir to run on: ex /srv/data/slxa/incoming/myrun123AAXX/Images\n$USAGE";
-d $runDir || die "DIR NOT FOUND. specify dir to run on: ex /srv/data/slxa/incoming/myrun/Images\n$USAGE";

#my @files = glob("$runDir/L00?/C*/*.tif");
my @files = glob("$runDir/L001/C14*/*.tif");
mkdir($workDir);
foreach my $image (@files)
{
	$semaphore->down();
	threads->new(\&downSample, $image);
}

for my $laneNumber (1..8)
{
	for my $cycleNumber (1..200)
	{
		$semaphore->down();
		threads->new(\&combine, $laneNumber, $cycleNumber);
	}
}

exit;

sub downSample
{
	my $tif = shift @_;
	$tif =~ /L00(\d)\/C(\d+)\.\d\/s_\d_(\d+)_(\w)\.tif/;
	my $lane = $1;
	my $cycle = $2;
	my $tile = $3;
	my $base = $4;

	$cycle = sprintf "%04d", $cycle;
	$tile = sprintf "%04d", $tile;
	my $outputJpeg =  "$workDir/$lane" . "_$cycle" . "_$tile" . "_$base" . ".jpg";
	if(!(-e $outputJpeg))
	{
		system("convert $tif -resize 10% $outputJpeg");
	}
	$semaphore->up();
}

sub combine
{
	my $lane = shift @_;
	my $cycle = shift @_;
	$cycle = sprintf "%04d", $cycle;
	my @cycleFiles = glob("$workDir/$lane" . "_$cycle" . "_*_a.jpg");
	if($#cycleFiles > 2)
	{
		my $cycleString = "";
		for my $i (0..59)
		{
			my $tile = sprintf "%04d", (60 - $i);
			my $nextTile =  sprintf "%04d", (61+$i);
			$cycleString .= "$workDir/$lane" . "_$cycle" . "_$tile" . "_a.jpg $workDir/$lane" . "_$cycle" . "_$nextTile" . "_a.jpg ";	
		}
		my $montExec = "montage $cycleString -tile 2x60 -geometry +0+0 $workDir/montage_$lane" . "_$cycle" . ".jpg";
		print "$montExec\n";
		system $montExec;
	}
	$semaphore->up();
}
