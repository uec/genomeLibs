#!/usr/bin/perl
@ARGV || die "usage: qseq2fastq.pl s_7_1_*qseq.txt\n";

for $qseq (@ARGV)
{
	open(IN,"<$qseq");
	while(<IN>)
	{
		@fields = split(/\t/);
		die "$qseq contains weird number of columns\n" if $#fields != 10;
		if($fields[10] > 0)
		{
			my $read = join("_", @fields[0..1]) . ":" . join(":", @fields[2..5]) . "#". join("/", @fields[6..7]);
			my $seq = $fields[8];
			$seq =~ s/\./N/g;
			my $qual = $fields[9];
			print "\@$read\n$seq\n\+$read\n$qual\n";
		}
	}
}

