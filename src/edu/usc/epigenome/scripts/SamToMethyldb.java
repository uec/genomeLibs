package edu.usc.epigenome.scripts;

import net.sf.samtools.*;

import java.io.File;
import java.sql.Connection;
import java.sql.DatabaseMetaData;
import java.sql.DriverManager;
import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;

import edu.usc.epigenome.genomeLibs.PicardUtils;


public class SamToMethyldb {

	final private static String TABLE = "methylCGsRich_tumor_chr11";
	// Any way to have ENUM wildcards?
	final private static String GET_BY_COORD_SQL_FW = "select * from " + TABLE + " WHERE chromPos = ? AND strand = '+';";
	final private static String GET_BY_COORD_SQL_REV = "select * from " + TABLE + " WHERE chromPos = ? AND strand = '-';";
	final private static String INSERT_BY_COORD_SQL_FW = "INSERT INTO " + TABLE + " VALUES (?,'+',?,?,?,?,?,?,?);";
	final private static String INSERT_BY_COORD_SQL_REV = "INSERT INTO " + TABLE + " VALUES (?,'-',?,?,?,?,?,?,?);";

	/**
	 * object vars
	 */
	protected Connection fConn;
	protected PreparedStatement fGetByCoordPrepFw;
	protected PreparedStatement fGetByCoordPrepRev;

	protected PreparedStatement fInsertByCoordPrepFw;
	protected PreparedStatement fInsertByCoordPrepRev;

	
	/**
	 * @param args
	 */
	@Option(name="-minConv",usage="minimum number of converted cytosines required")
	protected int minConv = 1;
//	@Option(name="-numCycles",usage="Number of cycles to track")
//	protected int numCycles = 100;
//	@Option(name="-outputReads",usage=" Outputs one line per read (default false)")
//	protected boolean outputReads = false;
	@Option(name="-useCpgsToFilter",usage=" Use CpGs and CpHs to filter if true, otherwise just CpHs (default false)")
	protected boolean useCpgsToFilter = false;
	@Option(name="-minMapQ",usage="minimum mapping quality (default 0)")
	protected int minMapQ = 0;
	@Option(name="-dryRun",usage=" Don't write output file")
	protected boolean dryRun = false;
	@Option(name="-debug",usage=" Debugging statements (default false)")
	protected boolean debug = false;

	
	// receives other command line parameters than options
	@Argument
	private List<String> stringArgs = new ArrayList<String>();

	
	
	/**
	 * @param args
	 */
	public static void main(String[] args)
	throws Exception
	{
		new SamToMethyldb().doMain(args);
	}

	public void doMain(String[] args)
	throws Exception {

		CmdLineParser parser = new CmdLineParser(this);
		// if you have a wider console, you could increase the value;
		// here 80 is also the default
		parser.setUsageWidth(80);
		try
		{
			parser.parseArgument(args);
			if (stringArgs.size() != 1) throw new CmdLineException("No input file specified");
		}
		catch (CmdLineException e)
		{
			System.err.println(e.getMessage());
			// print the list of available options
			parser.printUsage(System.err);
			System.err.println();
			return;
		}
		
		
		File inputSamOrBamFile = new File(stringArgs.get(0));

		if (!this.dryRun) this.setupDb();

		
		final SAMFileReader inputSam = new SAMFileReader(inputSamOrBamFile);
		inputSam.setValidationStringency(SAMFileReader.ValidationStringency.SILENT);
		int recCounter = 0;
		int usedCounter = 0;
		int filteredOutCounter = 0;
		record: for (final SAMRecord samRecord : inputSam) {
			
			// Filter low qual
			int mapQual = samRecord.getMappingQuality();
			boolean unmapped = samRecord.getReadUnmappedFlag();
			if (unmapped || (mapQual < minMapQ))
			{
				continue record;
			}
			
			String seq = PicardUtils.getReadString(samRecord, true);

			recCounter++;
			if ((recCounter % 1E5)==0) System.err.println("On new record #" + recCounter); 

			
			try
			{
				String ref = PicardUtils.refStr(samRecord, true);

				if (seq.length() != ref.length())
				{
					System.err.println("SeqLen(" + seq.length() + ") != RefLen(" + ref.length() + ")");
					System.err.println(seq + "\n" + ref);
				}
				//System.err.println(seq + "\n" + ref);

				boolean negStrand = samRecord.getReadNegativeStrandFlag();
				int	onRefCoord = (negStrand) ? samRecord.getUnclippedEnd() : samRecord.getAlignmentStart();
					
					
				int numConverted = 0;
				int convStart = Integer.MAX_VALUE;
				int seqLen = Math.min(seq.length(), ref.length());
				for (int i = 0; i < seqLen; i++)
				{
					char refi = ref.charAt(i);
					char seqi = seq.charAt(i);
					
					if ((seqi != '0') && PicardUtils.isCytosine(i,ref))
					{
						boolean cpg = PicardUtils.isCpg(i,ref);
						boolean conv = PicardUtils.isConverted(i,ref,seq);

						if (conv && (this.useCpgsToFilter || !cpg)) numConverted++;

						// If this is the first legal one , note it
						if ((convStart==Integer.MAX_VALUE) && (numConverted>=this.minConv) )
						{
							convStart = i;
						}

						
						if (cpg)
						{
							if (i<convStart)
							{
								// In the non-conversion filter zone
								filteredOutCounter++;
								//System.err.printf("Rec %d\tpos=%d\n",recCounter,i);
							}
							else
							{
								// Past the non-conversion filter, use it
								usedCounter++;
							}
							
							if (!this.dryRun) this.incrementDbCpg(onRefCoord, negStrand, seqi, i<convStart);

						}
						

						
					} // IsCytosine

					boolean oppositeCpg = PicardUtils.isOppositeCpg(i,ref);
					if (oppositeCpg)
					{
						if (!this.dryRun) this.incrementDbOppositeCpg(onRefCoord, negStrand, seqi);
					}

					
					// Increment genomic coord position
					if (refi == '-')
					{
						// It's a deletion in reference, don't advance
					}
					else
					{
						int inc = (negStrand) ? -1 : 1;
						onRefCoord += inc;
					}
					
				} // i (pos within read)
			}
			catch (Exception e)
			{
				System.err.println("-----------------------------------------");
				System.err.println("Couldn't handle seq #" + recCounter);
				System.err.println(seq);
				e.printStackTrace(System.err);
				System.err.println("-----------------------------------------");

				// Revert DB
				if (!this.dryRun)
				{
					fConn.rollback();
					fConn.close();
				}

				System.exit(1);
				

			}

		} // record

		
		double frac = (double)filteredOutCounter/((double)usedCounter+(double)filteredOutCounter);
		System.err.printf("Lost %f%% due to non-converion filter\n%d CpGs filtered for non-conversion, %d CpGs used (MinConv=%d,UseCpgs=%s)\n",
				frac*100.0, filteredOutCounter, usedCounter, this.minConv, String.valueOf(this.useCpgsToFilter));
		System.err.printf("Found %d reads total\n", recCounter);
		
		inputSam.close();
		
		if (!this.dryRun) this.cleanupDb();

		
	}

	
	protected void incrementDbCpg(int chromPos, boolean negStrand, char seqChar, boolean nonconvFilter) 
	throws Exception
	{
		int totalReads = 0, cReads = 0, tReads = 0, cReadsNonconvFilt = 0, agReads = 0;
		
		switch (seqChar)
		{
		case 'N':
		case '0':
			break;
		case 'A':
		case 'G':
			agReads = 1;
			totalReads = 1;
			break;
		case 'T':
			tReads = 1;
			totalReads = 1;
			break;
		case 'C':
			if (nonconvFilter) cReadsNonconvFilt = 1; else cReads = 1;
			totalReads = 1;
			break;
		default:
			throw new Exception("Can't recognize seq char: " + seqChar);
		}

		incrementDb(chromPos, negStrand, totalReads, cReads, cReadsNonconvFilt, tReads, agReads, 0 , 0);
	}
	
	protected void incrementDbOppositeCpg(int chromPos, boolean negStrand, char seqChar) 
	throws Exception
	{
		int aReadsOpposite = 0, totalReadsOpposite = 0;
		
		// Strand is reversed since we want relative to cytosine
		negStrand = !negStrand;
		
		switch (seqChar)
		{
		case 'N':
		case '0':
			break;
		case 'A':
			aReadsOpposite = 1;
			totalReadsOpposite = 1;
			break;
		case 'G':
		case 'T':
		case 'C':
			totalReadsOpposite = 1;
			break;
		default:
			throw new Exception("Can't recognize seq char: " + seqChar);
		}

		incrementDb(chromPos, negStrand, 0, 0, 0, 0, 0, aReadsOpposite , totalReadsOpposite);
	}
	
	public static void incrementShort(ResultSet rs, int col, short incVal)
	throws SQLException
	{
		short val = rs.getShort(col);
		rs.updateShort(col, (short)(val + incVal));
	}
	
	protected void incrementDb(int chromPos, boolean negStrand, int totalReads, int cReads,
			int cReadsNonconversionFilt, int tReads, int agReads, int aReadsOpposite, int totalReadsOpposite)
	throws SQLException
	{
		System.out.printf("pos=%d\tstrand=%c\ttotalReads=%d\tcReads=%d\tcReadsNonconv=%d\ttReads=%d\tagReads=%d\taReadsOpp=%d\ttotalReadsOpp=%d\n",
				chromPos,  negStrand ? '-' : '+',  totalReads,  cReads,
				 cReadsNonconversionFilt,  tReads,  agReads,  aReadsOpposite,  totalReadsOpposite);
		
		
		// First see if we already have this one
		PreparedStatement prep = (negStrand) ? this.fGetByCoordPrepRev : this.fGetByCoordPrepFw;
		prep.setInt(1, chromPos);
		ResultSet rs = prep.executeQuery();
		
		if (rs.next())
		{
			// Already in DB
			System.out.printf("%d%c: Already in database\n", chromPos, negStrand?'-':'+');
			incrementShort(rs, 3, (short)totalReads);
			incrementShort(rs, 4, (short)cReads);
			incrementShort(rs, 5, (short)cReadsNonconversionFilt);
			incrementShort(rs, 6, (short)tReads);
			incrementShort(rs, 7, (short)agReads);
			incrementShort(rs, 8, (short)totalReadsOpposite);
			incrementShort(rs, 9, (short)aReadsOpposite);
			
			rs.updateRow();

		}
		else
		{
			// First time seen
			System.out.printf("%d%c: first time encountered\n", chromPos, negStrand?'-':'+');
			prep = (negStrand) ? this.fInsertByCoordPrepRev : this.fInsertByCoordPrepFw;
			prep.setInt(1, chromPos);
			prep.setShort(2, (short)totalReads);
			prep.setShort(3, (short)cReads);
			prep.setShort(4, (short)cReadsNonconversionFilt);
			prep.setShort(5, (short)tReads);
			prep.setShort(6, (short)agReads);
			prep.setShort(7, (short)totalReadsOpposite);
			prep.setShort(8, (short)aReadsOpposite);

			prep.executeUpdate();
		}

		
	}
	
	protected void setupDb()
	throws Exception
	{
		String connStr = "jdbc:mysql://localhost/cr?user=benb";
		Class.forName("com.mysql.jdbc.Driver").newInstance();
		System.err.println("Getting connection for " + connStr);
		fConn = DriverManager.getConnection(connStr);
		
		fConn.setAutoCommit(false);

		// Make them type scroll sensitive so we can increment on the fly without autocommit
		fGetByCoordPrepFw = fConn.prepareStatement(GET_BY_COORD_SQL_FW, ResultSet.TYPE_SCROLL_SENSITIVE, ResultSet.CONCUR_UPDATABLE);
		fGetByCoordPrepRev = fConn.prepareStatement(GET_BY_COORD_SQL_REV, ResultSet.TYPE_SCROLL_SENSITIVE, ResultSet.CONCUR_UPDATABLE);

		fInsertByCoordPrepFw = fConn.prepareStatement(INSERT_BY_COORD_SQL_FW);
		fInsertByCoordPrepRev = fConn.prepareStatement(INSERT_BY_COORD_SQL_REV);
	
		
//		try {
//			DatabaseMetaData dmd = fConn.getMetaData();
//			if (dmd.supportsResultSetConcurrency(
//					ResultSet.TYPE_FORWARD_ONLY, ResultSet.CONCUR_UPDATABLE)) {
//				System.err.println("Updatable sets are supported.");
//			} else {
//				// Updatable result sets are not supported
//				System.err.println("Updatable sets are NOT supported.");
//			}
//
//			if (dmd.supportsTransactions()) {
//				System.err.println("Transactions are supported.");
//			} else {
//				System.err.println("Transactions are NOT supported.");
//			}
//
//		   
//		   } catch (SQLException e) {
//		    }
		
	}
	
	protected void cleanupDb()
	throws Exception
	{
		fConn.commit();
		fConn.close();
	}

	
//
//	static int boolToInd(boolean bool)
//	{
//		return ( bool ? 0 : 1 );
//	}

	
	
}
