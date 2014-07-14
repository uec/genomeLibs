#!/usr/bin/perl

#GLOBALS
$findpeaks = "/home/uec-00/shared/production/software/VancouverShortR/fp4/FindPeaks.jar";
$JAVA = "/home/uec-00/shared/production/software/java/default/bin/java";

#INPUTS
$inputFile = $ARGV[0] || die "specifiy input BAM file";
$inputFile =~ /^(.+)\.bam$/ || die "specficy BAM file";
$prefix = $1;
$dist = $ARGV[1] || die "specify freg length (ex: 200)";
my $pid = fork();

#$cmd = "java -jar $findpeaks -no_warning -name $prefix -input $inputFile -output ./ -aligner sam -dist_type 0 $dist -landerwaterman 0.001 -duplicatefilter  -auto_threshold 0.001";
#system($cmd);
#exit;

if($pid == 0)
{
	#$cmd = "$JAVA -jar $findpeaks -eff_frac 0.7 -no_warning -name $prefix -input $inputFile -output ./ -aligner sam -dist_type 1 $dist -iterations 5 -duplicatefilter -auto_threshold 0.01";
	#print "$cmd\n";
	#system("$cmd");
}
else
{
	$cmd = "$JAVA -jar $findpeaks -no_warning -name $prefix\_raw -input $inputFile -output ./ -aligner sam -dist_type 1 $dist -duplicatefilter";
	print "$cmd\n";
	system("$cmd");
	waitpid($pid,0);
}

