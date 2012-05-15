#!/usr/bin/perl -w
# Mon Jul 12 14:03:06 PDT 2010

my $usage = "wrap_novo.pl refFa.fa read1.fq [read2.fq] [outfile name]";
die "$usage\n" unless (@ARGV >= 3);

my $novo = "/home/uec-00/shared/production/software/novoalign/default/novoalign";
my $SAMTOOLS = "/home/uec-00/shared/production/software/samtools/samtools";
my $PICARD = "/home/uec-00/shared/production/software/picard/default";
my $JAVA = "/home/uec-00/shared/production/software/java/default/bin/java";

$refNovo = $ARGV[0];
#catch fasta passed instead of novo index.
$refNovo =~ s/\.fa$/\.novoindex\.bis/;

#and then set the fasta for later use with picard
$refFa = $refNovo;
$refFa =~ s/novo.+$/fa/;

$read1 = $ARGV[1];
$read2 = $ARGV[2] if(@ARGV == 4);
$outfile = $ARGV[$#ARGV];

die "$outfile already exists" if (-e $outfile);

$outfileSAM = $outfile . ".sam";

# check existence of ref genome
die "reference does not exist.\n" if (! -e $refNovo);

# check existence of read sequence files
die "need read sequence files\n" unless ( -e $read1 );

#aln end 1
my $cmd = "$novo -d $refNovo -f $read1 $read2 -F ILMFQ -b 4 -o SAM > $outfileSAM";
runcmd($cmd);

runcmd("$JAVA -Xmx14g -jar $PICARD/SortSam.jar VALIDATION_STRINGENCY=SILENT SORT_ORDER=coordinate INPUT=$outfileSAM OUTPUT=$outfile\.sorted.bam");
runcmd("$JAVA -Xmx14g -jar $PICARD/ReorderSam.jar VALIDATION_STRINGENCY=SILENT REFERENCE=$refFa INPUT=$outfile\.sorted.bam OUTPUT=$outfile");

sub runcmd
{
        my $cmd = shift @_;
        print STDERR "$cmd\n";
        system($cmd);
}
