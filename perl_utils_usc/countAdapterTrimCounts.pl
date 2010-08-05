#!/usr/bin/perl

$output = shift || die "Specifiy output file";
die "$output already exists" if -e $output;

for $file (@ARGV)
{
	if(-e $file)
	{
		open(IN, "<$file");
		while(<IN>)
		{
			chomp;
			$counts[$_]++;
		}	
	}
}

open(OUT,">$output") || die "output error";
for $i (0..$#counts)
{
	my $count = $counts[$i] || 0;
	print OUT "adapterTrimCount $i $count\n";
}
