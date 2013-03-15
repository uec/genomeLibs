#!/usr/bin/perl


my $JAVA = "/home/uec-00/shared/production/software/java/default/bin/java -Xmx1995m -classpath /home/uec-00/shared/production/software/genomeLibs/apps-live.jar:/home/uec-00/shared/production/software/genomeLibs/biojava-live.jar:/home/uec-00/shared/production/software/genomeLibs/bytecode.jar:/home/uec-00/shared/production/software/genomeLibs/commons-cli.jar:/home/uec-00/shared/production/software/genomeLibs/commons-collections-2.1.jar:/home/uec-00/shared/production/software/genomeLibs/commons-dbcp-1.1.jar:/home/uec-00/shared/production/software/genomeLibs/commons-math-1.1.jar:/home/uec-00/shared/production/software/genomeLibs/commons-pool-1.1.jar:/home/uec-00/shared/production/software/genomeLibs/demos-live.jar:/home/uec-00/shared/production/software/genomeLibs/genomeLibs.jar:/home/uec-00/shared/production/software/genomeLibs/jgrapht-jdk1.5.jar:/home/uec-00/shared/production/software/genomeLibs/junit-4.4.jar:/home/uec-00/shared/production/software/genomeLibs/charts4j-1.2.jar:/home/uec-00/shared/production/software/genomeLibs/heatMap.jar:/home/uec-00/shared/production/software/genomeLibs/UscKeck.jar";
my $fqfa = "/home/uec-00/shared/production/software/perl_utils_usc/fastqToFasta.pl";
my $output =  shift @ARGV || die "file not found";
my $input = join(" ", @ARGV);


for $i (3,5,10)
{
	my $nmers = $i . "mers";
	$output =~ s/\d+mers/$nmers/;
	my $cmd = "cat $input | $fqfa | $JAVA edu.usc.epigenome.scripts.FastaToNmerCounts -nmer $i | head -n 1000 > $output";  
	print STDERR "$cmd\n";
	system($cmd);
}

