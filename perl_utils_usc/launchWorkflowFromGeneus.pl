#!/usr/bin/perl
use LWP::UserAgent;
use strict;
our @files;
my $processURL = $ARGV[0];
$ARGV[0] =~ /http:\/\// || die "usage: program.pl http://processURL";

#get flowcell ID from remote xml using given URI
my $browser = LWP::UserAgent->new;
$browser->credentials(
  'epilims.usc.edu:8080',
  'GLSSecurity',
  'zack' => 'genzack'
);
my $response = $browser->get($processURL);
my $processContent =  $response->content;
$processContent =~ /Flowcell S\/N\"\>([\d\w]+)\</ || die "could not fine flowcell serial in remote XML";
my $flowcell = $1;

print "$flowcell\n";

#find the sequence files on our storage servers
&findFiles("/storage/gastorage2/slxa/incoming");
&findFiles("/storage/polycomb/slxa/incoming");
exit() if scalar @files < 1;

#make remote project dir:
my $mkdirTries = 1;
my $remoteDir;
my $mkdirOutput;
#do
#{
	$remoteDir = "projects/$flowcell/run$mkdirTries";
	$mkdirOutput = `ssh -l ramjan hpc-uec.usc.edu mkdir $remoteDir 2>&1`;
	$mkdirTries++;
#} while {$mkdirOutput =~ /File exists/i && $mkdirTries < 5};

die "error creating project dir at hpcc, already exists?" if length($mkdirOutput) > 10;
print $mkdirOutput;

#copy inputs to remote project dir
for my $f (@files)
{
	 `scp $f ramjan\@hpc-uec.usc.edu:$remoteDir`;
}

#launch job
`ssh -l ramjan hpc-uec \"cd $remoteDir; java -cp /home/uec-00/shared/production/software/ECWorkflow/ECWorkFlow.jar:/home/uec-00/shared/production/software/ECWorkflow/pegasus.jar edu.usc.epigenome.workflow.AlignPileupWorkflow -pbs $processURL\"`;


#find all seq result files in a given server mount point
sub findFiles
{
	my $startDir = shift @_;
	chdir($startDir);
	my $findCmd = `find *$flowcell* -size +100M -name  \"*s_?_sequence.txt\" 2>/dev/null`;
	#print $findCmd;
	my @hits = split(/\s*\n\s*/, $findCmd);

	#make sure we only have at least the 8 lanes, fail if more
	$#hits < 8 || die "aborting, too many input sequences found: $findCmd";
	for my $x (@hits)
	{
		if(length($x) > 10)
		{
			push @files, "$startDir/$x";
		}
	}
	#make sure we dont have same flowcell on multiple servers
	$#files < 8 || die "aborting, found same flowcell on multiple servers. ambiguous";
}
