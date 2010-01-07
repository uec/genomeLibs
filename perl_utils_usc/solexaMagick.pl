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

#synchronize
$semaphore->down(8);
$semaphore->up(8);

for my $laneNumber (1..8)
{
	for my $cycleNumber (1..200)
	{
		$semaphore->down();
		threads->new(\&combineCycle, $laneNumber, $cycleNumber);
	}
}

#synchronize
$semaphore->down(8);
$semaphore->up(8);


for my $laneNumber (1..8)
{
	$semaphore->down();
	&combineLane($laneNumber)	
}

#synchronize
$semaphore->down(8);
$semaphore->up(8);

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
		system("convert $tif -resize 20% $outputJpeg");
	}
	$semaphore->up();
}

sub combineCycle
{
	my $lane = shift @_;
	my $cycle = shift @_;
	$cycle = sprintf "%04d", $cycle;
	&combineCycleByBase($lane,$cycle,"a");
	&combineCycleByBase($lane,$cycle,"c");
	&combineCycleByBase($lane,$cycle,"g");
	&combineCycleByBase($lane,$cycle,"t");
	my @cycleBaseFiles = glob("$workDir/montage_$lane" . "_$cycle" . "_?.jpg");
	my $cycleBaseString = join(' ', sort(@cycleBaseFiles));
	my $montExec = "montage $cycleBaseString -tile 1x4 -geometry +0+0 $workDir/montage_$lane" . "_$cycle.jpg";
	print "$montExec\n";
	system $montExec;
	$semaphore->up();
}

sub combineCycleByBase
{
	my $lane = shift @_;
	my $cycle = shift @_;
	my $base = shift @_;
	$cycle = sprintf "%04d", $cycle;
	my @cycleFiles = glob("$workDir/$lane" . "_$cycle" . "_*_$base" . ".jpg");
	if($#cycleFiles > 2)
	{
		my $cycleString = "";
		for my $i (0..59)
		{
			my $tile = sprintf "%04d", (60 - $i);
			my $nextTile =  sprintf "%04d", (61+$i);
			$cycleString .= "$workDir/$lane" . "_$cycle" . "_$tile" . "_$base" . ".jpg $workDir/$lane" . "_$cycle" . "_$nextTile" . "_$base" . ".jpg ";	
		}
		my $montExec = "montage $cycleString -tile 2x60 -geometry +0+0 $workDir/montage_$lane" . "_$cycle" . "_$base" . ".jpg";
		print "$montExec\n";
		system $montExec;
	}
}

sub combineLane
{
		my $lane = shift @_;
		my @cycleFiles = glob("$workDir/montage_$lane" . "_????" . ".jpg");
		my $cycleString = join(' ', sort(@cycleFiles));
		my $numCycles = scalar @cycleFiles;
		my $montExec = "montage $cycleString -tile $numCycles" . "x1 -geometry +0+0 $workDir/montage_$lane" . ".jpg";
		print "$montExec\n";
		system $montExec;
		$semaphore->up();	
}

