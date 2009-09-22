#!/usr/bin/perl
use File::Spec;

$dir = $ARGV[0] || die "specify a dir to compress";
$cpus = $ARGV[1] || die "specify number of cpus to use";

-d $dir || die "dir not found";
$cpus =~ /^\d+$/ || die "cpus must be a number";

@files = split(/\n/,`find $dir -name \"*bz2\"`);

$tmpBase = $dir;
$tmpBase =~ /^.+\/([\d\w_\-]+?)$/;
$tmpBase = $1;

#split data into many pbs jobs
foreach $i  (0..$#files)
{
	$cmd[$i % $cpus] .= "bzip2 -d" . File::Spec->rel2abs($files[$i]) ."\n";
}


#create pbs files
foreach $i (0..$#cmd)
{
	my $header = "#PBS -q laird\n";
	$header .= "#PBS -l walltime=48:00:00\n";
	my $outfileName = "pbsDecompress" . $tmpBase . ".$i.sh";
	mkdir("pbsDecompress_" . $tmpBase);
	chdir("pbsDecompress_" . $tmpBase);
	open(OUT, ">$outfileName");
	print OUT $header;
	print OUT $cmd[$i];
	close OUT;
	system("qsub $outfileName");
	chdir("..");

}


