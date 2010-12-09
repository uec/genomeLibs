#!/usr/bin/perl
#cat s_7_1_*_qseq.txt | awk -F '\t' '{gsub(/\./,"N", $9); if ($11 > 0) print "@"$1"_"$2":"$3":"$4":"$5":"$6"#"$7"/"$8"\n"$9"\n""+"$1$2$3$4$5$6$8"\n"$1 0}'

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

