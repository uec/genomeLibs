#!/usr/bin/env perl

use strict;
use File::Basename;
use Getopt::Long;
use File::Temp qw/ tempfile /;

my $USAGE = "removeRedundantEntries.pl in1.gtf in2.gtf";


my $centeredSize = 0;
GetOptions ('centeredSize=i'=>\$centeredSize)|| die "$USAGE\n";

die "$USAGE\n" unless (@ARGV>0);

foreach my $fn (@ARGV)
{
	die "Can't read $fn\n" unless (open(FIN,$fn)); 


    my ($name, $path, $suf) = fileparse($fn, qr/\.g[tf]f/); #qr/\.[^.]*/);
	my $outFn = "${name}.nodups.gtf";

	my $fieldsByKey = {};
	
	while (my $line = <FIN>)
	{
		chomp $line;
		my @f = split(/\t/,$line);
		my $gene = $f[8];
		#print STDERR "Line has ". scalar(@f) . " fields\n";
		

		my $key = toKey(\@f);
		#print STDERR "key=$key\n";
		my $oldFields = $fieldsByKey->{$key};
		
		my $update = 0;
		if (!$oldFields)
		{
			$update = 1;
		}
		else
		{
			my $oldGene = @$oldFields[8];
		
			if (!$oldGene && $gene)
			{
				$update = 1;			
			} 		
		}
		
		if ($update)
		{
			$fieldsByKey->{$key} = \@f;	
		}
	}
	close (FIN); 


	# Output
	die "Can't write to $outFn\n" unless (open(FOUT, ">$outFn"));
	foreach my $key (sort(keys(%$fieldsByKey)))
	{
		my $f = $fieldsByKey->{$key};
		print FOUT join("\t",@$f)."\n";
	}
	
	close(FOUT);	
}

sub toKey
{
	my ($f) = @_;
	return join("_",@$f[0],@$f[3],@$f[4],@$f[6]);
}
