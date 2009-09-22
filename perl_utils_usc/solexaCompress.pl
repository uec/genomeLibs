#!/usr/bin/perl
use File::Spec;

$dir = $ARGV[0] || die "specify a dir to compress";
$cpus = $ARGV[1] || die "specify number of cpus to use";

-d $dir || die "dir not found";
$cpus =~ /^\d+$/ || die "cpus must be a number";

@files = split(/\n/,`find $dir -size +512k`);
@files = grep(!/htm.*$/i,@files);
@files = grep(!/png$/i,@files);
@files = grep(!/gz$/i,@files);
@files = grep(!/zip$/i,@files);
@files = grep(!/jpg$/i,@files);
@files = grep(!/bzip2$/i,@files);
@files = grep(!/bz.*$/i,@files);

$tmpBase = $dir;
$tmpBase =~ /^.+\/([\d\w_\-]+?)$/;
$tmpBase = $1;

#split data into many pbs jobs
foreach $i  (0..$#files)
{
	$cmd[$i % $cpus] .= "bzip2 " . File::Spec->rel2abs($files[$i]) ."\n";
}


#create pbs files
foreach $i (0..$#cmd)
{
	my $header = "#PBS -q laird\n";
	$header .= "#PBS -l walltime=48:00:00\n";
	my $outfileName = "pbsCompress" . $tmpBase . ".$i.sh";
	mkdir("pbsCompress_" . $tmpBase);
	chdir("pbsCompress_" . $tmpBase);
	open(OUT, ">$outfileName");
	print OUT $header;
	print OUT $cmd[$i];
	close OUT;
	system("qsub $outfileName");
	chdir("..");

}


