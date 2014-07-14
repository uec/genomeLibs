#!/usr/bin/perl
use File::Basename;
die "must be a dir" unless -d $ARGV[0];
chdir($ARGV[0]);

$mask = "-name \"s_*sequence.[0-9]*.*\" ! -name \"*tophat*\"";



system("find $mask | xargs du -ch");
#echo "find -name \"*sequence\.[0-9]*\.*\" | xargs rm"
system("find $mask -exec rm {} \\;"); 


@mdupBams= glob("*.mdups.bam");
for my $mdup (@mdupBams)
{
	if($mdup =~ /\.mdups\.bam$/)
	{
		my $rawbam = basename($mdup);
		$rawbam =~ s/\.mdups\.bam$/\.bam/;
		if(-s $rawbam && -s $mdup && $rawbam ne $mdup && length($mdup) > length($rawbam))
		{
			
			print ("rm $rawbam\n");
			print ("rm $rawbam\.bai\n");
		}
	}

}

