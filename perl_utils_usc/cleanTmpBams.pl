#!/usr/bin/perl
use strict;
use File::Basename;

my $dir = $ARGV[0] || die "what dir to look in?";
-d $dir || die "must be a dir";
open FOUND, "find $dir -mtime +3 -mtime -60 -name \"*.bam\" |";



while(my $tbam = <FOUND>)
{
        chomp $tbam;
        my $dirname  = dirname($tbam);
        my @bams = glob("$dirname/*.bam");
	my @toDelete = @bams;
	my @mdups = grep(/mdups/,@bams); 
	@mdups = sort { -s $b <=> -s $a } @mdups;
	@toDelete = grep(!/mdups/,@bams);
	@toDelete = grep(!/unmapped/,@toDelete);
	@toDelete = grep(!/pbs/,@toDelete);
	@toDelete = sort { -s $b <=> -s $a } @toDelete;

	if(-e $mdups[0] && -e $toDelete[0] && $mdups[0] ne $toDelete[0] && -s $mdups[0] >= -s $toDelete[0] )
	{
		print "KEEPING  $_\n" for (@mdups);
		print "DELETING $_\n" for (@toDelete);
		system("rm $_") for (@toDelete);
	}
	print "\n";

}
