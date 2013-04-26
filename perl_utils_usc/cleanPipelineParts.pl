#!/usr/bin/perl
die "must be a dir" unless -d $ARGV[0];
chdir($ARGV[0]);

$mask = "\"*.[0-9]*.fastq\" -o -name \"s_*sequence.[0-9]*.*\" ! -name \"*tophat*\"";



system("find $mask | xargs du -ch");
#echo "find -name \"*sequence\.[0-9]*\.*\" | xargs rm"
system("find $mask -exec rm {} \\;"); 

@mdups = glob("*mdups.bam");

for my $bam (@mdups)
{
	$reg = $bam;
	$reg =~ s/\.mdups//;

	print ("rm $reg") if (-e $reg && -e $bam && $reg !~ /mdups/);


}

