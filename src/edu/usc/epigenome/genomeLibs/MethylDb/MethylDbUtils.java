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

import org.biojava.bio.seq.StrandedFeature;

import edu.usc.epigenome.genomeLibs.ChromScores.ChromScoresArray;
import edu.usc.epigenome.genomeLibs.ChromScores.ChromScoresArrayInt;
import edu.usc.epigenome.genomeLibs.ChromScores.ChromScoresFast;
import edu.usc.epigenome.genomeLibs.FeatDb.FeatDbQuerier;

public class MethylDbUtils {
	public static final List<String> CHROMS =
		Arrays.asList(
//				"chrY", "chrM",
				"chr1","chr2","chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", 
        		"chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", 
        		"chr19", "chr20", "chr21", "chr22", "chrX");
	
	public static final List<String> SMALL_CHROMS =
		Arrays.asList(
//				"chr6", "chr7", "chr8","chr9",
        		"chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", 
        		"chr19", "chr20", "chr21", "chr22");

	public static final List<String> TEST_CHROMS =
	Arrays.asList("chr11");


	
	public static final List<String> SUMMARY_FEATURES1 = 
		Arrays.asList("IMR90_PMDs","TJ_Prmtrs_oriented","K27me3_Ku2008_Prmtrs_oriented","K27me3_Ku2008");
	public static final List<String> SUMMARY_FEATURES2 = 
		Arrays.asList("LINE","exon","IMR90_PMDs","TJ_Prmtrs_oriented","K27me3_Ku2008_Prmtrs_oriented","K27me3_Ku2008");
	
	protected static Connection cConn = null; 
	protected static Map<String,PreparedStatement> cPreps = new HashMap<String,PreparedStatement>();


	public static double fetchMeanExpression(String chr, String refseqId, String sqlExpression)
	throws Exception
	{
		
		boolean refseqFormat = ((refseqId != null) && refseqId.matches("^[A-Z][A-Z]_.*"));
		//System.err.printf("RefseqId=%s, matches format=%s\n", refseqId, ""+refseqFormat);
		if (!refseqFormat) return Double.NaN;
		
		//System.err.println("RefseqId=" + refseqId);
		
		setupDb();
		
		String table = "infiniumExpr_" + chr;
		String sql = String.format("SELECT %s FROM %s exp WHERE exp.refseqId = ?;", sqlExpression, table);
		PreparedStatement prep = MethylDbUtils.getPrep(sql);
		prep.setString(1, refseqId);
		ResultSet rs = prep.executeQuery();
		
		//System.err.println("Fetching expr for : " + refseqId);
		int numFound = 0;
		double out = 0.0;
		while (rs.next())
		{
			//System.err.println("\tGot expression val: " + rs.getDouble(1));
			out += rs.getDouble(1);
			numFound++;
		}
		
		//System.err.println("numFound = " + numFound);
		if (numFound == 0)
		{
			out = Double.NaN;
		}
		else if (numFound>1)
		{
			Logger.getLogger(Logger.GLOBAL_LOGGER_NAME).log(Level.FINE,
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

	/**
	 * @param rec
	 * @param consensusCpg
	 * @return return null if not enough cpgs
	 */
	public static String bedLine(String chr, int s, int e, String strand, double meth)
	{
		

		int score;
		String color = methToColor(meth);
		if (Double.isNaN(meth) || Double.isInfinite(meth))
		{
			score = 0;
		}
		else
		{
			score = (int)Math.round(100.0 * meth);
		}
//		else if (meth > 0.5)
//		{
//			score = (int)Math.round(100.0 * ((meth-0.5)*2.0));
//		}
//		else
//		{
//			score = (int)Math.round(100.0 * ((0.5-meth)*2.0));
//		}

		String out = String.format("%s\t%d\t%d\t%s\t%d\t%s\t%d\t%d\t%s\t", 
				chr,
				s,
				e,
				String.format("r%d%s", s,strand),
				score,
				strand,  // strand information gets in the way of display
				s,
				e,
				color
				);
		return out;
	}

	public static String methToColor(double meth)
	{
		String color="255,255,204"; // 

		if (Double.isNaN(meth))
		{
		}
		else if (meth > 0.5)
		{
			int dec = (int)Math.floor(10.0 * ((meth-0.5)*2.0));
			// Red
//			switch (dec)
//			{
//			case 0: color = "49,49,49"; break;
//			case 1: color = "55,44,44"; break;
//			case 2: color = "61,38,38"; break;
//			case 3: color = "66,33,33"; break;
//			case 4: color = "72,27,27"; break;
//			case 5: color = "78,22,22"; break;
//			case 6: color = "83,16,16"; break;
//			case 7: color = "89,11,11"; break;
//			case 8: color = "94,5,5"; break;
//			case 9: case 10: color = "100,0,0"; break;
//			default: System.err.println("Got illegal color decile: " + dec); System.exit(1); break;
//			}
			switch (dec)
			{
			case 0: color = "153,153,153"; break;
			case 1: case 2: case 3: color = "255,153,153"; break;
			case 4: case 5: case 6: case 7: color = "255,102,102"; break;
			case 8: case 9: case 10: color = "255,0,0"; break;
			default: System.err.println("Got illegal color decile: " + dec); System.exit(1); break;
			}
		}
		else
		{
			int dec = (int)Math.floor(10.0 * ((0.5-meth)*2.0));
			// Green
//			switch (dec)
//			{
//			case 0: color = "49,49,49"; break;
//			case 1: color = "44,55,44"; break;
//			case 2: color = "38,61,38"; break;
//			case 3: color = "33,66,33"; break;
//			case 4: color = "27,72,27"; break;
//			case 5: color = "22,78,22"; break;
//			case 6: color = "16,83,16"; break;
//			case 7: color = "11,89,11"; break;
//			case 8: color = "5,94,5"; break;
//			case 9: case 10: color = "0,100,0"; break;
//			default: System.err.println("Got illegal color decile: " + dec); System.exit(1); break;
//			}
			switch (dec)
			{
			case 0: color = "153,153,153"; break;
			case 1: case 2: case 3: color = "153,255,153"; break;
			case 4: case 5: case 6: case 7: color = "102,255,102"; break;
			case 8: case 9: case 10: color = "0,255,0"; break;
			default: System.err.println("Got illegal color decile: " + dec); System.exit(1); break;
			}
		}
	
		return color;
	}
	
	public static ChromScoresFast[] chromScoresReadCounts(MethylDbQuerier params, String chr, String tablePrefix, String inGenome)
	throws Exception
	{
		return chromScoresReadCounts( params,  chr,  tablePrefix,  inGenome, 0, (int)2.8E8);
	}
	
	public static ChromScoresFast[] chromScoresReadCounts(MethylDbQuerier params, String chr, String tablePrefix, String inGenome, int chromS, int chromE)
	throws Exception
	{
		
		// Setup  array
		ChromScoresFast[] scores = new ChromScoresFast[2];
		scores[0] = new ChromScoresArrayInt(inGenome); // FW STRAND
		scores[1] = new ChromScoresArrayInt(inGenome); // REV STRAND
		
		int MINCOORD = chromS;
		int MAXCOORD = (chromE>0) ? chromE : (int)2.8E8;
		int STEP = (int)1E6;
		
		int numSeen = 0;
		for (int s = MINCOORD; s < MAXCOORD; s += STEP)
		{
			// Get iterator
			CpgIteratorMultisample it = chromCpgIteratorMultisample(params, chr, tablePrefix, s, s+STEP-1);

			// Populate array
			while (it.hasNext())
			{
				if ((numSeen%1E5)==0) System.err.printf("On CpG #%d\n",numSeen);
				try
				{
					Cpg[] cpgs = it.next();
					Cpg cpg = cpgs[0];
					double count = cpg.totalReads;
					ChromScoresFast strandScores = (cpg.getStrand() == StrandedFeature.NEGATIVE) ? scores[1] : scores[0];
					//			strandScores.addScore(chr, cpg.chromPos, count);
					strandScores.addScore(chr, cpg.chromPos, 1);
					numSeen++;
				}
				catch (Exception e)
				{
					System.err.println("Skipping CpG: ");
					e.printStackTrace();
				}
			}
		}
		
		return scores;
	}

	public static ChromScoresFast[] chromScoresMethLevels(MethylDbQuerier params, String chr, String tablePrefix, String inGenome)
	throws Exception
	{
		return chromScoresMethLevels( params,  chr,  tablePrefix,  inGenome, 0, (int)2.8E8);
	}

	public static ChromScoresFast[] chromScoresMethLevels(MethylDbQuerier params, String chr, String tablePrefix, String inGenome, int chromS, int chromE)
	throws Exception
	{
		// Setup  array
		ChromScoresFast[] scores = new ChromScoresFast[2];
		scores[0] = new ChromScoresArrayInt(inGenome);  // CpG positions [0,1]
		scores[1] = new ChromScoresArrayInt(inGenome);  // Meth totals
		
		int MINCOORD = chromS;
		int MAXCOORD = (chromE>0) ? chromE : (int)2.8E8;
		int STEP = (int)1E6;
		
		int numSeen = 0;
		for (int s = MINCOORD; s < MAXCOORD; s += STEP)
		{

			// Get iterator
			CpgIteratorMultisample it = chromCpgIteratorMultisample(params, chr, tablePrefix, s, s+STEP-1);

			// Populate array
			int lastPos = -1;
			while (it.hasNext())
			{
				if ((numSeen%1E5)==0) System.err.printf("On CpG #%d\n",numSeen);
				try
				{
					Cpg[] cpgs = it.next();
					Cpg cpg = cpgs[0];
					int pos = cpg.chromPos;
					if (pos == lastPos) System.err.printf("Why did we see coord %d twice?\n",pos);

					double meth = cpg.fracMeth(params.useNonconversionFilter);
					if (meth>1.0) System.err.printf("Meth>1.0 (%.3f)\n",meth);

					// If it's a minus strand one, slide it back to the + coord.
					//if (cpg.getStrand() == StrandedFeature.NEGATIVE) pos--;

					scores[0].addScore(chr, pos, 1.0);
					scores[1].addScore(chr, pos, 100.0 * meth);
					lastPos = pos;
					numSeen++;
				}
				catch (Exception e)
				{
					System.err.println("Skipping CpG: ");
					e.printStackTrace();
				}
			}
		}
		
		return scores;
	}

	private static CpgIteratorMultisample chromCpgIteratorMultisample(MethylDbQuerier params, String chr, 
			String tablePrefix, int rangeStart, int rangeEnd)
	throws Exception
	{
		params.clearRangeFilters();
		params.addRangeFilter(chr, rangeStart, rangeEnd);
//		params.addRangeFilter(chr,(int)14E6,(int)15E6);
		
		CpgIteratorMultisample it = null;
		it = new CpgIteratorMultisample(params, Arrays.asList(tablePrefix));
		
		params.clearRangeFilters();
		return it;
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
