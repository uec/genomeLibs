#!/usr/bin/perl


=AUTHOR
Dennis T. Maglinte, dmaglinte@gmail.com
Copyright (c) 2011 Epicenter Software


=DESCRIPTION
Parse data from *.fastq.tophat_hits.bam.cufflinks_transcripts.gtf produced by Zack's RNASeq pipeline. File to be readable by Genetrix. Parse data
to temp array of array. After temp AoA is created, sort by chromosome and then sort by chr start position. Input file is not truly sorted
although it appears to be sorted alphabetically.
 
note: For GTF to BED conversion, subtract 1 from all starts


=INPUT
chr1	Cufflinks	transcript	757758	759190	1000	-	.	gene_id "CUFF.921803"; transcript_id "CUFF.921803.1"; FPKM "0.5730178534"; frac "0.707513"; conf_lo "0.000000"; conf_hi "2.105276"; cov "1.558399";
chr1	Cufflinks	exon	757758	759190	1000	-	.	gene_id "CUFF.921803"; transcript_id "CUFF.921803.1"; exon_number "1"; FPKM "0.5730178534"; frac "0.707513"; conf_lo "0.000000"; conf_hi "2.105276"; cov "1.558399";
chr1	Cufflinks	transcript	757995	759190	592	-	.	gene_id "CUFF.921803"; transcript_id "CUFF.921803.2"; FPKM "0.3396204828"; frac "0.292487"; conf_lo "0.000000"; conf_hi "1.553317"; cov "0.923643";
chr1	Cufflinks	exon	757995	758018	592	-	.	gene_id "CUFF.921803"; transcript_id "CUFF.921803.2"; exon_number "1"; FPKM "0.3396204828"; frac "0.292487"; conf_lo "0.000000"; conf_hi "1.553317"; cov "0.923643";
chr1	Cufflinks	exon	758194	759190	592	-	.	gene_id "CUFF.921803"; transcript_id "CUFF.921803.2"; exon_number "2"; FPKM "0.3396204828"; frac "0.292487"; conf_lo "0.000000"; conf_hi "1.553317"; cov "0.923643";


=OUTPUT
example:
1	656298	6571926	7.0	+	2|0.15|0.50:0:2934:1826:73:3489:406:2344:1073|0.50:0:2934:1826:73:3489:100:20


layout:
chr	start	stop	total_abundance	strand	conf_lo	conf_hi	n isoforms|p|isoform(relative abundance|5' offset|exon #1 length|intron #1 length|exon #2 length| intron #2...)


=cut


# Global
$dir_in = "";
$dir_out = "";
my $file_in = "HSB_107_AMY_L.1.nocontam.fastq.tophat_hits.bam.cufflinks_transcripts";
my $fastq = $dir_in . $file_in . ".gtf";
my $txt_out = $dir_out . $file_in . ".genetrix.txt";
my $num_chr = 24;
my $x = $num_chr-1;
#length = #stop-start +1


#====================================================================================================
# Clear temp file
#====================================================================================================
open (OUT, ">$txt_out")|| die("Cannot open file \"$txt_out\"\n\n");
print OUT "track type=isoform chr_col=1 start_col=2 end_col=3 val_col=4 strand_col=5 label_col=6 label_col=7 label_col=8 gapvalue=missing\n";


#====================================================================================================
# Read in all splice variants from fastq file and parse them into a temp file
#====================================================================================================
#print "Reading in file: $fastq...\n";

# none
open (FASTQ, " < $fastq") || die("Cannot open file \"$fastq\"\n\n");

my @temp_gene; # Store info at gene level. To be reset after printed.
my @temp_isoforms; # Store isoforms for given gene. To be reset after printed.
my $bool_gene = 0; # Keep track that gene info was already stored.
my $bool_iso = 0;
my $prev_end = 0;
my $temp_iso = "";
my $prev_CUFF = ""; # Keep track of what gene we are on. Same CUFF ID for all isoforms of a given gene.

my @AoA_genes;
my @num_genesInchr; # Number of genes in each chromosome.

# Default number of genes in each chr is 0
for (my $i=0; $i<$num_chr; $i++) {
	$num_genesInchr[$i] = 0;
}

#print "Parsing file.\n";
while(<FASTQ>) {
	my ($seq, $start, $end, $total_abundance, $strand, $transcript_status, $transcript_num, $conf_lo, $conf_hi, $cov, $rel_abundance);
	my @line = split(" ",$_);
	my $cur_CUFF = substr($line[9],6, length($line[9])-8);
	
	$seq = $line[0];
	if ($seq eq "chrM") {
		goto nextline;
	}
	
	# Get numeric version of chromsome number, minus "chr"
	$seq = substr($seq, 3);
	if ($seq =~ /^[\dXY]+$/) {
		$seq=~ s/X/$x/;	
		$seq=~ s/Y/$num_chr/;	
	}
	$start = $line[3]-1; #BED
	$end = $line[4];
	$total_abundance = substr($line[13],1,length($line[13])-3);
	$strand = $line[6];
	$conf_lo = substr($line[17],1,length($line[17])-3);
	$conf_hi = substr($line[19],1,length($line[19])-3);
	$cov = substr($line[21],1,length($line[21])-3);
	$transcript_status = $line[2];
	$transcript_num = substr($line[11], length($line[9])-1, length($line[11])-(length($line[9])-1)-2);
	
	
	# Moved on to new gene so delete temp arrays.
	# Don't print transcript info until CUFF ID changes. Unique CUFF ID per gene.	
	# Print out gene info (applies to all isoforms of a given gene)
	if ($prev_CUFF ne $cur_CUFF && $prev_CUFF ne "") {
		if ($bool_iso) {
			push (@temp_isoforms, $temp_iso);
		}

		my $iso = "";
		for (my $i=0; $i<@temp_isoforms; $i++) {
			$iso .= "$temp_isoforms[$i]|"
		}
		chop $iso;
		push (@temp_gene,$iso); # Save isoform info as a single line in the array 
		push @AoA_genes, [@temp_gene]; # Store gene in array of array for sorting
		$num_genesInchr[$temp_gene[0]-1]++;
		
	
		# Reset variables for new gene 
		undef @temp_gene;
		undef @temp_isoforms;
		$temp_iso = "";
		$bool_iso = $bool_gene = $prev_end = 0;
	}

	# Check if we are looking at a new transcript/isoform
	if ($transcript_status eq "transcript") {
		
		if ($bool_iso) {
			push (@temp_isoforms, $temp_iso);
			
			# Isoform saved in temp array so reset string
			$temp_iso = "";
			$bool_iso = $prev_end = 0;
		}
		
		if (!$bool_gene) {
			# Haven't started looking at this gene yet
			$temp_gene[0] = $seq;
			$temp_gene[1] = $start;
			$temp_gene[2] = $end;
			$temp_gene[3] = $total_abundance;
			$temp_gene[4] = $strand;
			$temp_gene[5] = $conf_lo;
			$temp_gene[6] = $conf_hi;
			$bool_gene = 1;
		} else {
			# Have looked at this gene already, but are looking at a different isoform
			if ($end > $temp_gene[2]) {
				# New end is bigger so replace old end.
				$temp_gene[2] = $end;
			}
			
			# TEMP*****************************************
			$temp_gene[3] += $total_abundance;
			#**********************************************
		}
	# Start reading in exon info
	} elsif ($transcript_status eq "exon") {
		$temp_isoforms[0] = $transcript_num; # number of isoforms
		
		if (!$bool_iso) {
			# First exon of a given isoform
			#TEMP**********************************
			$temp_isoforms[1] = 0;#$cov; # p-value
			#**************************************
			$rel_abundance = substr($line[17],1,length($line[17])-3);
		
			my $offset5prime = $start - $temp_gene[1];
			my $exon_length = $end - $start;
			
			$temp_iso = $rel_abundance . ":" . $offset5prime . ":" . $exon_length;
			$bool_iso = 1;
			$prev_end = $end;
		} else {
			# Next or last exon of a given isoform
			my $intron_length = $start - $prev_end;
			my $exon_length = $end - $start;
			
			$temp_iso .= ":" . $intron_length . ":" . $exon_length;
			$prev_end = $end;
		}	
	}
	$prev_CUFF = $cur_CUFF;
	nextline:
}

if ($bool_iso) {
	push (@temp_isoforms, $temp_iso);
}

my $iso = "";
for (my $i=0; $i<@temp_isoforms; $i++) {
	$iso .= "$temp_isoforms[$i]|"
}
chop $iso;
push (@temp_gene,$iso); # Save isoform info as a single line in the array 
push @AoA_genes, [@temp_gene]; # Store gene in array of array for sorting
$num_genesInchr[$temp_gene[0]-1]++;


close (FASTQ);



#====================================================================================================
# Sort array of array of genes
#====================================================================================================

#print "Sorting file.\n";

# Sort first on chromosome (first element in array)
@AoA_genes = sort { $a->[0] <=> $b->[0] } @AoA_genes;

my $c = 0;
for ($i=0; $i<$num_chr; $i++) {
	
	# Get all genes for one chromosome into a temp AoA
	my @AoA_temp;
	for (my $j=$c; $j<$c+$num_genesInchr[$i]; $j++) {
		push @AoA_temp, [@{$AoA_genes[$j]}];
	}
	$c += $num_genesInchr[$i];

	# Sort temp AoA by chromosome start position (second element in array)
	@AoA_temp = sort { $a->[1] <=> $b->[1] } @AoA_temp;

	
#====================================================================================================
# Print array of array of genes
#====================================================================================================
	
	for (my $j=0; $j <$#{AoA_temp}+1; $j++) {
		# Substitute back X and Y for chromosome if needed
		$AoA_temp[$j][0] =~ s/$x/X/;
		$AoA_temp[$j][0] =~ s/$num_chr/Y/;
		
		print OUT "$AoA_temp[$j][0]\t$AoA_temp[$j][1]\t$AoA_temp[$j][2]\t$AoA_temp[$j][3]\t$AoA_temp[$j][4]\t$AoA_temp[$j][5]\t$AoA_temp[$j][6]\t$AoA_temp[$j][7]\n";
	}	
	undef @AoA_temp;
}

close (OUT);



#print "DONE\n\n";