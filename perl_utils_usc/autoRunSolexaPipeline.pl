#!/usr/bin/perl

use LWP::Simple;
use XML::Simple;
use Data::Dumper;
use File::Basename;
use Cwd;
use strict;

my $SERVER = "gastorage1";

my $content = get("http://www.epigenome.usc.edu/gareports/flowcell.php?xml") || die "Couldn't get it!";
my $xml = XMLin($content, KeyAttr=>["serial","lane"]);
my @readyDirs;

chdir("/srv/$SERVER/slxa/incoming");
foreach my $file (glob("*/Basecalling_Netcopy_complete.txt"))
{
        my $dir = dirname($file);
        my @processedDirs = glob("$dir/Data/Intensities/BaseCalls/GERALD* $dir/MAKELOG*");
        if(!$processedDirs[0])
        {
                push @readyDirs, $dir;
        }
}

foreach my $dir (@readyDirs)
{
        &doPipeline($dir);
}

sub doPipeline
{
        print "\n\n";
        use vars qw($xml);
        my $dir = shift @_;
        $dir =~ /(\w\w\w\w\wAAXX)/;

        my $curDir = getcwd();

        my %organismList;
        my $flowcell = $1 || die "flowcell name not in dir name";
        foreach my $i (1..8)
        {
                $organismList{$xml->{flowcell}->{"$flowcell"}->{sample}->{$i}->{organism}} .= $i;
        }

        my %analysisList = %organismList;
        my %genomeList = %organismList;

        my $config = "#auto generated config for flowcell $flowcell\n";
        $config .= "EMAIL_LIST ramjan\@usc.com benbfly\@gmail.com\n";
        $config .= "WEB_DIR_ROOT http://gastorage2.usc.edu/slxa_runs/incoming/\n";
        $config .= "12345678:SRF_ARCHIVE_REQUIRED yes\n";
        $config .= "12345678:SRF_QCAL yes\n";

        for my $org (keys %analysisList)
        {
                my $analysis = ($org =~  /^phi/i) ? "ANALYSIS eland_extended" : "ANALYSIS sequence";
                $config .= "$analysisList{$org}:$analysis\n"
        }

        $config .= "ELAND_MULTIPLE_INSTANCES 8\n";
        for my $org (keys %genomeList)
        {
                my $genome;
                if($org =~  /^phi/i) { $genome = "/srv/data/slxa/GENOMES/phi-X174/" ; }
                elsif($org =~ /^Homo/i) { $genome = "/srv/data/slxa/GENOMES/hg18_unmasked/"; }
                else { print "Unknown Organism or geneus data not found, skipping $dir\n";return; }
                $config .= "$analysisList{$org}:ELAND_GENOME $genome\n"
        }

        #write config
        print $config;
        my $configFileName = "config.sequence.auto.txt";
        chdir($dir);
        open(my $configFile, ">$configFileName") || die "FAILED, couldnt write to config File\n";
        print $configFile $config;

        #create MakeFile
        my $geraldCmd = "/opt/GAPipeline-1.6.0a9/bin/GERALD.pl $configFileName --ok_to_use_legacy_gerald  --EXPT_DIR ./Data/Intensities/BaseCalls --make >\& MAKELOG1";
        print "running:\n$geraldCmd\n";
        system($geraldCmd);

        #run MakeFile
        my @geraldDir = glob("Data/Intensities/BaseCalls/GER*");
        if(-d $geraldDir[0])
        {
                my $makeCmd = "nohup make -j 8 recursive&";
                chdir($geraldDir[0]);
                print("running:\n$makeCmd\n");
                system($makeCmd);
        }
        else
        {
                print "couldn't locate GERALD Makefile";
        }
        chdir($curDir);
}