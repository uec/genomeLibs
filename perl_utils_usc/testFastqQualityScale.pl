#!/usr/bin/perl
$file = $ARGV[0] || die;

checkFile($ARGV[0],$_) for (10,100,10000,100000,1000000);
print "sanger phred33\n";



sub checkFile
{
	my $file = shift @_;
	my $checks = shift @_;
	my $size = -s $file || die;

	open(IN,$file);

	for (my $i=0; $i < $size; $i += ($size/$checks))
	{
		seek(IN,int($i),0);
		my @lines;
		for my $j (1..8)
		{
			my $line = <IN>;
			push @lines, $line;
		}
		
		my $j=0;
		$j++ while (!($lines[$j] =~ /^\@/ && $lines[$j+2] =~ /^\+/) && $j < 9);
		if($lines[$j+3] =~ /[KLMNOPQRSTUVWXYZ\[\\\]\^\_\`abcdefgh]+/)
		{
			print STDERR int($i), $lines[$j+3];
			print "illumina phred64\n";
			exit;
		}
		if($lines[$j+3] =~ /[\!\"\#\$\%\&\'\(\)\*\+\,\-\.\/0123456789\:\;\<\=\>\?]+/)
		{
			print STDERR int($i), $lines[$j+3];
			print "sanger phred33\n";
			exit;
		}
	}
	close(IN);
}
