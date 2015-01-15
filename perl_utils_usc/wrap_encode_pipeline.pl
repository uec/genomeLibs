#!/usr/bin/perl
use strict;
use File::Basename;
use Getopt::Long;
use threads;
use threads::shared;
use Data::Dumper;

use lib dirname (__FILE__);
use EpigenomeUtils;

#the following should be set according to your system
my $BEDTOOLS_PATH = "$SOFTWAREROOT/bedtools/default/bin";
my $R = "$SOFTWAREROOT/R/release/R-3-x86_64/bin/R";
my $SPP = "$SOFTWAREROOT/SPP/default/run_spp_nodups.R";
my $IDR = "$SOFTWAREROOT/SPP/idrCode/batch-consistency-analysis.r";

#since we have a single R install shared amongst many people, I've installed the spp libraries locally, so I set R_LIBS
$ENV{R_LIBS} = "$SOFTWAREROOT/SPP/default/Rlibs";
$ENV{PATH} .= ":" .  dirname($SAMTOOLS);

my $MAXMEM_KB_PER_JOB = 6291456;
my $MAX_CONCURRANT_JOBS = `cat /proc/meminfo | grep MemTotal`;
$MAX_CONCURRANT_JOBS =~ s/\D//g;
$MAX_CONCURRANT_JOBS = int($MAX_CONCURRANT_JOBS / $MAXMEM_KB_PER_JOB);
print STDERR "MAX_CONCURRANT_JOBS : $MAX_CONCURRANT_JOBS\n";

my $step :shared = 1;
my $activeJobs :shared = 0;

my (@treatedBams,@controlBams,$verbose);
GetOptions ("treated=s" => \@treatedBams, "control=s" => \@controlBams, "verbose" => \$verbose);
die "wrap_encode_pipeline.pl -control control.bam -treated sample1.bam -treated sample2.bam" unless @treatedBams && @controlBams;

#check existance of all files

map {die "$_ NOT FOUND\n" unless -e $_} @treatedBams;
map {die "$_ NOT FOUND\n" unless -e $_} @controlBams;

my %cleanNames;
$cleanNames{cleanName($treatedBams[$_])} = "S" . ($_ + 1) for (0..$#treatedBams);
$cleanNames{cleanName($controlBams[$_])} = "C" .  ($_ + 1) for (0..$#controlBams);

open (MAPPING,">step_000.Sample.Mapping");
print MAPPING "labelling $_ as $cleanNames{$_}\n" for sort keys %cleanNames;
close MAPPING;

#create TagAlign Files
my (@controlTags, @treatedTags);
my (@controlTagThreads, @treatedTagThreads);
push @controlTagThreads, threads->create(\&toTagAlign,$_) for @controlBams;
push @treatedTagThreads, threads->create(\&toTagAlign,$_) for @treatedBams;

push @controlTags, $_->join() for @controlTagThreads;
push @treatedTags, $_->join() for @treatedTagThreads;

#pool controls
my $pooledControl = poolTA(@controlTags);

#pool treated
my $pooledSample = poolTA(@treatedTags);


#create Psuedoreps
my (@psuedoReps1, @psuedoReps2, @pseudoThreads);
push @pseudoThreads, threads->create(\&splitFile,$_) for @treatedTags;

for my $t (@pseudoThreads)
{
	my ($r1,$r2) = $t->join();
	push @psuedoReps1, $r1;
	push @psuedoReps2, $r2;
}

#pool pseudoReps
my $pooledPseudoRep1 = poolTA(@psuedoReps1);
my $pooledPseudoRep2 = poolTA(@psuedoReps2);


#call replicates

#thread lists
my (@treatedRePPeaksThreads,@treatedPooledPeaksThreads,@pseudoPeaksThreadsR1,@pseudoPeaksThreadsR2,@pseudoPeaksThreadsPool1,@pseudoPeaksThreadsPool2);

#call reps
push @treatedRePPeaksThreads, threads->create(\&spp,$_,$pooledControl) for @treatedTags;

#call pooled 
push @treatedPooledPeaksThreads, threads->create(\&spp,$pooledSample,$pooledControl);

#call psuedo reps
push @pseudoPeaksThreadsR1, threads->create(\&spp,$_,$pooledControl) for @psuedoReps1;
push @pseudoPeaksThreadsR2, threads->create(\&spp,$_,$pooledControl) for @psuedoReps2;
push @pseudoPeaksThreadsPool1, threads->create(\&spp,$pooledPseudoRep1,$pooledControl);
push @pseudoPeaksThreadsPool2, threads->create(\&spp,$pooledPseudoRep2,$pooledControl);


#wait for SPP threads to finish
my (@treatedRePPeaks,@treatedPooledPeaks,@pseudoPeaksR1,@pseudoPeaksR2,@pseudoPeaksPool1,@pseudoPeaksPool2);
push @treatedRePPeaks, $_->join() for @treatedRePPeaksThreads;
push @treatedPooledPeaks, $_->join() for @treatedPooledPeaksThreads;
push @pseudoPeaksR1, $_->join() for @pseudoPeaksThreadsR1;
push @pseudoPeaksR2, $_->join() for @pseudoPeaksThreadsR2;
push @pseudoPeaksPool1, $_->join() for @pseudoPeaksThreadsPool1;
push @pseudoPeaksPool2, $_->join() for @pseudoPeaksThreadsPool2;
@pseudoPeaksR1 = sort @pseudoPeaksR1;
@pseudoPeaksR2 = sort @pseudoPeaksR2;

#idr on reps
my @idrThreads;
my %seenSampleRePeaks;
for my $peak1 (@treatedRePPeaks)
{
	$seenSampleRePeaks{"$peak1\-$peak1"} = 1;
	for my $peak2 (@treatedRePPeaks)
	{
		push @idrThreads,threads->create(\&idr,$peak1,$peak2) unless $seenSampleRePeaks{"$peak1\-$peak2"} || $seenSampleRePeaks{"$peak2\-$peak1"};
		$seenSampleRePeaks{"$peak1\-$peak2"} = 1;
		$seenSampleRePeaks{"$peak2\-$peak1"} = 1;
	}
}


#idr on pseudo pools
push @idrThreads,threads->create(\&idr,$pseudoPeaksPool1[0],$pseudoPeaksPool2[0]);

#idr on indiv pseudo pools
push @idrThreads,threads->create(\&idr,$pseudoPeaksR1[$_],$pseudoPeaksR2[$_]) for (0..$#pseudoPeaksR1);

#Wait for IDR threads to finish
$_->join() for @idrThreads;


##########################################################################

sub toTagAlign
{
	my $input = shift @_;
	return $input if $input =~ /lign.gz$/i;

	my $workdir = "step_" . getStep() . ".toTagAlign_" . cleanName($input); 
	die "Output already exists, are you running this on top of old results, delete the old stuff first (ex: rm -r step_* )\n" if -e $workdir;
	mkdir $workdir;
	my $output = cleanName($input) . ".tagAlign.gz";
	$output = basename($output);
	
	runcmd("$SAMTOOLS view -b -F 1548 -q 30 $input | $BEDTOOLS_PATH/bamToBed -i stdin | awk 'BEGIN{FS=\"\\t\";OFS=\"\\t\"}{\$4=\"N\"; print \$0}' | gzip -c > $workdir/$output") if $input =~ /\.bam$/i;
	runcmd("cat $input | gzip -c > $workdir/$output") if $input =~ /\.bed$/ || $input =~ /\.tagAlign$/i ;
	return "$workdir/$output";
}

sub poolTA
{
	my @input = @_;
	my $inputList = join(" ",@input);
	print STDERR "pooling { $inputList }\n";
	my @output = map {basename $_} @input;
	my $outputFile = cleanName(join("_",@output));
	my $workdir = "step_" . getStep() . ".poolSamples_" . $outputFile; 
	mkdir $workdir;
	$outputFile .= ".joined.tagAlign.gz";
	runcmd("zcat $inputList | gzip -c > $workdir/$outputFile");
	return "$workdir/$outputFile";
}

sub spp
{
	my $input = shift @_;
	my $control = shift @_;
	my $label = cleanName($input) . "_VS_" . cleanName($control);
	print STDERR "Running SPP ON $input $control\n";
	$label =~ s/\///g;
	my $workdir = "step_" . getStep() . ".SPP_callpeaks_$label"; 
	mkdir $workdir;
	runcmd("$R" . "script $SPP -c=$input -i=$control -npeak=300000 -odir=$workdir -savr -savp -rf -out=$workdir/phantomPeakStatsReps.tab");
	runcmd("gzip -d $workdir/*Peak.gz");
	return glob("$workdir/*Peak");
}

sub splitFile
{
	my $input = shift @_;
	my $lines = `zcat $input  | wc -l`;
	chomp $lines;
	$lines = int(($lines + 1) / 2);
	my $output = cleanName($input) . ".pr";
	my $workdir = "step_" . getStep() . ".split_pseudoreplicate_" . cleanName($input); 
	mkdir $workdir;
	runcmd("zcat $input | shuf | split -d -l $lines - $workdir/$output"); 
	my @files = glob("$workdir/*pr*");
	runcmd("mv $_ $_\.tagAlign") for @files;
	runcmd("gzip $_\.tagAlign") for @files;
	return glob("$workdir/*tagAlign.gz");
}

sub idr
{
	my $input = shift @_;
	my $control = shift @_;
	print STDERR "performing IDR analysis on $input $control\n";
	my $output = cleanName($input) . "__IDR-VS__" .  cleanName($control);
	my $workdir = "step_" . getStep() . ".IDR_analysis_" . $output; 
	mkdir $workdir;
	runcmd("$R" . "script $IDR $input $control -1 $workdir/$output 0 F signal.value");
	return glob("$workdir/*");
}


sub runcmd
{
	getResources();
	my $cmd=shift @_;
	my $caller=(caller(1))[3];
	print STDERR "$caller\t$cmd\n";
	$verbose ? system($cmd) : `$cmd`;
	releaseResources();
}

sub cleanName
{
	my $input = shift @_;
	my $output = basename($input);
	$output =~ s/gz//g;
	$output =~ s/\.bam//g;
	$output =~ s/\.bed//g;
	$output =~ s/tagAlign//g;
	$output =~ s/\.\.+/\./g;
	$output =~ s/\._/_/g;
	$output =~ s/\.$//;
	$output =~ s/_$//;
	for my $key (%cleanNames)
	{
		#print STDERR "relabbeling $key to $cleanNames{$key} \n " if $output eq $key && $cleanNames{$key};
		return $cleanNames{$key} if $output eq $key && $cleanNames{$key};

	}
	return $output;
}

sub getStep
{
	my $i = $step++;
	return "00" . $i if $i < 10;
	return "0" . $i if $i < 100; 
	return $i;
}

#since thread pools mod was not installed, and dont have root, this func will only return when resources are avail.
sub getResources
{
	sleep(rand(10)) while($activeJobs >= $MAX_CONCURRANT_JOBS);

	{
		lock($activeJobs); 
		$activeJobs++;
	}
}
sub releaseResources
{
	{
		lock($activeJobs); 
		$activeJobs--;
	}
}
