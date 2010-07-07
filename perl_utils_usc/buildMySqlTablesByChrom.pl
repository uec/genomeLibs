#!/usr/bin/env perl

use strict;
use File::Basename;
use Getopt::Long;
use File::Temp qw/ tempfile /;

my $USAGE = "buildMySqlTablesByChrom.pl [--noclobberTable|--templateFile schema.sql] [--pipeConverterScript pipe.pl] [--db cr] [--replaceStr TEMPLATE_CHR] table1_chr1.txt table1_chr1.txt ...";

my $delim = "\t";
my $templateFn = 0;
my $replaceStr = "TEMPLATE_CHR";
my $pipe = 0;
my $pipePrefix = 0;
my $clobberTable = 1;
my $db = "cr";
GetOptions ('clobberTable!'=>\$clobberTable, 'pipeConverterScript=s'=>\$pipe, 'pipePrefix=s'=>\$pipePrefix, 'templateFile=s' => \$templateFn, 'replaceStr=s' => \$replaceStr) || die "$USAGE\n";

print STDERR "clobber=$clobberTable\targv=" . scalar(@ARGV) . "\n";
die "$USAGE\n" unless ((!$clobberTable || $templateFn) && (@ARGV>0));

foreach my $tabfn (@ARGV)
{
    my ($name, $path, $suf) = fileparse($tabfn, qr/\.[^.]*/);

 	if ($clobberTable)
 	{
		die "Can't read $templateFn\n" unless (open(TF,$templateFn));
		my @templateLines = ();
		while (my $line = <TF>)
		{
			push(@templateLines,$line);
		}
		close(TF);
	
		my $tmpfn = makeTableFile(\@templateLines, $replaceStr, $name);
		
		runCmd("cat ${tmpfn} | mysql ${db}");
		unlink($tmpfn);
		#runCmd("echo \"drop table ${name};\" | mysql ${db}"); #----REMOVE
	}

	# Do we need to pipe the file through a converter?  Since the file has to have the same name as the 
	# table, we have to move it.
	if ($pipe)
	{
		$name =~ /(chr\w+)/;
		my $nameChr = $1;
		my $tabfnNew = ($pipePrefix) ? "${pipePrefix}_${nameChr}.txt" : "${name}.converted.txt";
		my $cmd = "${pipe} < $tabfn > $tabfnNew";
		runCmd($cmd);
		$tabfn = $tabfnNew;
	}


	my $deleteSec = ($clobberTable) ? " --delete " : "";	
	my $cmd = "mysqlimport ${deleteSec} --use-threads=2 --lock-tables --replace --force --local --verbose ${db} ${tabfn}";
	runCmd($cmd);
}

sub runCmd
{
	my ($cmd, $testOnly) = @_;
	print STDERR $cmd."\n";
	print STDERR `$cmd`."\n" unless ($testOnly);
}

sub makeTableFile
{
	
	my ($templateLines, $replaceStr, $withStr) = @_;
	
	my ($fh, $filename) = tempfile("buildMySqlTablesByChrom.pl.XXXXXX");
	
	for (my $i = 0; $i < scalar(@$templateLines); $i++)
	{
		my $l = @{$templateLines}[$i];
		$l =~ s/$replaceStr/$withStr/g;
		print $fh $l;
	}
	
	close($fh);
	return $filename;
}