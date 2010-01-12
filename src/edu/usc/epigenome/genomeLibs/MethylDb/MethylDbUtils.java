package edu.usc.epigenome.genomeLibs.MethylDb;

import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;

import edu.usc.epigenome.genomeLibs.FeatDb.FeatDbQuerier;

public class MethylDbUtils {
	public static final List<String> CHROMS =
		Arrays.asList("chr1","chr2","chr4", "chr5", "chr6", "chr7", "chr8", "chr9", 
        		"chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", 
        		"chr19", "chr20", "chr21", "chr22", "chrX"); // "chr3");//, "chrY", "chrM");//  Why is chrom 3 missing??
	
	public static final List<String> TEST_CHROMS =
	Arrays.asList("chr11","chr12");


	
	public static final List<String> SUMMARY_FEATURES1 = 
		Arrays.asList("IMR90_PMDs","TJ_Prmtrs_oriented","K27me3_Ku2008_Prmtrs_oriented","K27me3_Ku2008");
	public static final List<String> SUMMARY_FEATURES2 = 
		Arrays.asList("LINE","exon","IMR90_PMDs","TJ_Prmtrs_oriented","K27me3_Ku2008_Prmtrs_oriented","K27me3_Ku2008");
	
	protected static Connection cConn = null; 
	protected static Map<String,PreparedStatement> cPreps = new HashMap<String,PreparedStatement>();


	public static double fetchMeanExpression(String chr, String refseqId, String sqlExpression)
	throws Exception
	{
		setupDb();
		
		String table = "infiniumExpr_" + chr;
		String sql = String.format("SELECT %s FROM %s exp WHERE exp.refseqId = ?;", sqlExpression, table);
		PreparedStatement prep = MethylDbUtils.getPrep(sql);
		prep.setString(1, refseqId);
		ResultSet rs = prep.executeQuery();
		
		int numFound = 0;
		double out = 0.0;
		while (rs.next())
		{
			out = rs.getDouble(1);
			numFound++;
		}
		
		//System.err.println("numFound = " + numFound);
		if (numFound == 0)
		{
			out = Double.NaN;
		}
		else if (numFound>1)
		{
			Logger.getLogger(Logger.GLOBAL_LOGGER_NAME).log(Level.SEVERE,
					String.format("Found multiple expr rows with refseq %s: %d rows\n", refseqId,numFound));
			out = out / (double)numFound;
		}
		
		return out;
	}
	
	protected static void setupDb()
	throws Exception
	{
		if (cConn == null)
		{
			String connStr = FeatDbQuerier.connStr;
			Class.forName("com.mysql.jdbc.Driver").newInstance();
			System.err.println("Getting connection for " + connStr);
			cConn = DriverManager.getConnection(connStr);
		}
		
	}
	
	protected static void cleanupDb()
	throws Exception
	{
		cConn.close();
	}

	protected static PreparedStatement getPrep(String sql)
	throws SQLException
	{
		PreparedStatement prep = cPreps.get(sql);
		if (prep==null)
		{
			prep = cConn.prepareStatement(sql);
			cPreps.put(sql, prep);
			Logger.getLogger(Logger.GLOBAL_LOGGER_NAME).log(Level.SEVERE, "Making prepared statement: " + sql );
		}
		return prep;
	}

	
//	//	echo "select count(*),featType from features_chr1 GROUP BY featType;" |mysql cr > featTypes.txt
//	61841   exon
//	16659   exon_normHigh
//	12843   exon_normLow
//	1730    exon_tumDown
//	15314   exon_tumHigh
//	12843   exon_tumLow
//	787     exon_tumUp
//	3228    FANTOM4_HCP
//	4128    FANTOM4_LCP
//	602     IMR90_NotPMDs
//	1206    IMR90_PMDs
//	1524    K27me3_Ku2008
//	416     K27me3_Ku2008_Prmtrs_oriented
//	228504  LINE
//	250     RING1B_Ku2008
//	323448  SINE
//	870     TJ_GG_nonPrmtrExons
//	1254    TJ_NonPrmtrNonExon
//	2140    TJ_Prmtrs_oriented
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
