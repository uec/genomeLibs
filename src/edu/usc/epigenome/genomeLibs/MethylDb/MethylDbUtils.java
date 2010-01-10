package edu.usc.epigenome.genomeLibs.MethylDb;

import java.util.Arrays;
import java.util.List;

public class MethylDbUtils {
	public static final List<String> CHROMS =
		Arrays.asList("chr1","chr2","chr4", "chr5", "chr6", "chr7", "chr8", "chr9", 
        		"chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", 
        		"chr19", "chr20", "chr21", "chr22"); // "chr3", "chrX", "chrY", "chrM");  Why is chrom 3 missing??
	
	public static final List<String> TEST_CHROMS =
	Arrays.asList("chr11","chr12", "chr13");



	
	public static final List<String> SUMMARY_FEATURES2 = 
		Arrays.asList("TJ_GG_nonPrmtrExons","TJ_NonPrmtrNonExon","TJ_Prmtrs_oriented");
	public static final List<String> SUMMARY_FEATURES3 = 
		Arrays.asList("exon_normHigh", "exon_normLow", "exon_tumHigh", "exon_tumLow", "exonTumUp", "exonTumDown");
	public static final List<String> SUMMARY_FEATURES4 = 
		Arrays.asList("LINE","exon","IMR90_PMDs","TJ_Prmtrs_oriented","K27me3_Ku2008_Prmtrs_oriented","K27me3_Ku2008");
	public static final List<String> SUMMARY_FEATURES5 = 
		Arrays.asList("tss_500bp_flank_normHigh","tss_500bp_flank_normLow", "tss_500bp_flank_tumHigh","tss_500bp_flank_tumLow","tss_500bp_flank_tumUp","tumDown");
	public static final List<String> SUMMARY_FEATURES1 = 
		Arrays.asList("tx_normHigh","tx_normLow", "tx_tumHigh","tx_tumLow");
	
	//	echo "select count(*),featType from features_chr1 GROUP BY featType;" |mysql cr > featTypes.txt
//	61841   exon
//	16659   exon_normHigh
//	12843   exon_normLow
//	1730    exon_tumDown
//	15314   exon_tumHigh
//	12843   exon_tumLow
//	787     exon_tumUp
//	1614    FANTOM4_HCP
//	2064    FANTOM4_LCP
//	603     IMR90_PMDs
//	762     K27me3_Ku2008
//	208     K27me3_Ku2008_Prmtrs_oriented
//	114252  LINE
//	125     RING1B_Ku2008
//	161724  SINE
//	435     TJ_GG_nonPrmtrExons
//	627     TJ_NonPrmtrNonExon
//	1070    TJ_Prmtrs_oriented
//	6633    tss
//	6633    tss_1kb_flank
//	1495    tss_1kb_flank_normHigh
//	1203    tss_1kb_flank_normLow
//	152     tss_1kb_flank_tumDown
//	1376    tss_1kb_flank_tumHigh
//	1203    tss_1kb_flank_tumLow
//	91      tss_1kb_flank_tumUp
//	6633    tss_2kb_flank
//	1495    tss_2kb_flank_normHigh
//	1203    tss_2kb_flank_normLow
//	152     tss_2kb_flank_tumDown
//	1376    tss_2kb_flank_tumHigh
//	1203    tss_2kb_flank_tumLow
//	91      tss_2kb_flank_tumUp
//	6633    tss_500bp_flank
//	1495    tss_500bp_flank_normHigh
//	1203    tss_500bp_flank_normLow
//	152     tss_500bp_flank_tumDown
//	1376    tss_500bp_flank_tumHigh
//	1203    tss_500bp_flank_tumLow
//	91      tss_500bp_flank_tumUp
//	1495    tss_normHigh
//	1203    tss_normLow
//	152     tss_tumDown
//	1376    tss_tumHigh
//	1203    tss_tumLow
//	91      tss_tumUp
//	6633    tx
//	1495    tx_normHigh
//	1203    tx_normLow
//	152     tx_tumDown
//	1376    tx_tumHigh
//	1203    tx_tumLow
//	91      tx_tumUp

	
}
