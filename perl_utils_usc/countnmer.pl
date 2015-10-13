#!/usr/bin/perl
use File::Basename;
use lib dirname (__FILE__);
use EpigenomeUtils;


my $JAVA = "$SOFTWAREROOT/java/default/bin/java -Xmx1995m -classpath $SOFTWAREROOT/genomeLibs/apps-live.jar:$SOFTWAREROOT/genomeLibs/biojava-live.jar:$SOFTWAREROOT/genomeLibs/bytecode.jar:$SOFTWAREROOT/genomeLibs/commons-cli.jar:$SOFTWAREROOT/genomeLibs/commons-collections-2.1.jar:$SOFTWAREROOT/genomeLibs/commons-dbcp-1.1.jar:$SOFTWAREROOT/genomeLibs/commons-math-1.1.jar:$SOFTWAREROOT/genomeLibs/commons-pool-1.1.jar:$SOFTWAREROOT/genomeLibs/demos-live.jar:$SOFTWAREROOT/genomeLibs/genomeLibs.jar:$SOFTWAREROOT/genomeLibs/jgrapht-jdk1.5.jar:$SOFTWAREROOT/genomeLibs/junit-4.4.jar:$SOFTWAREROOT/genomeLibs/charts4j-1.2.jar:$SOFTWAREROOT/genomeLibs/heatMap.jar:$SOFTWAREROOT/genomeLibs/UscKeck.jar";
my $fqfa = "$SOFTWAREROOT/perl_utils_usc/fastqToFasta.pl";
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

