package edu.usc.epigenome.genomeLibs.MethylDb;

import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.PreparedStatement;
import java.sql.ResultSet;

import java.sql.SQLException;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;
import java.util.logging.Level;
import java.util.logging.Logger;


/**
 * @author benb
 *
 */
public class CpgIterator implements Iterator<Cpg> {

	// Class vars
	protected static Connection cConn = null; 
	protected static Map<String,PreparedStatement> cPreps = new HashMap<String,PreparedStatement>();
	
	// Object vars
	protected MethylDbQuerier params = null;
	protected ResultSet curRS = null;
	protected int curNumRows = -1;
	

	/**
	 * 
	 */
	public CpgIterator(MethylDbQuerier inParams)
	throws Exception 
	{
		super();
		this.init(inParams);
	}



	
	public int getCurNumRows() {
		return curNumRows;
	}




	public CpgIterator(String chrom, int startCoord, int endCoord, String inTablePrefix)
	throws Exception 
	{
		super();
		params = new MethylDbQuerier();
		if (inTablePrefix!=null) params.methylTablePrefix = inTablePrefix;
		params.addRangeFilter(chrom, startCoord, endCoord);
		this.init(params);
	}
		
	
	public CpgIterator(String chrom, String inTablePrefix)
	throws Exception 
	{
		super();
		params = new MethylDbQuerier();
		if (inTablePrefix!=null) params.methylTablePrefix = inTablePrefix;
		params.addRangeFilter(chrom);
		this.init(params);

	}

	public CpgIterator(MethylDbQuerier inParams, String inTablePrefix)
	throws Exception 
	{
		params = inParams;
		if (inTablePrefix!=null) params.methylTablePrefix = inTablePrefix;
		this.init(params);
	}	
	
	public int init(MethylDbQuerier inParams)
	throws Exception
	{
		this.params = inParams;

		// Check if we've started DB connection
		if (cConn == null) setupDb();
		PreparedStatement prep = CpgIterator.getPrep(inParams);
		CpgIterator.fillPrep(params, prep);
		Logger.getLogger(Logger.GLOBAL_LOGGER_NAME).fine("Starting query execute");
		curRS = prep.executeQuery();		
		curRS.last();
		int numRows = curRS.getRow();
		curRS.beforeFirst();
		Logger.getLogger(Logger.GLOBAL_LOGGER_NAME).fine("Finished query execute");
		
		this.curNumRows = numRows;
		return numRows;
	}
	
	
	
	public CpgIterator()
	throws Exception 
	{
//		throw new Exception("Whole genome Cpg iterator not yet supported");
//		super();
//		// Check if we've started DB connection
//		if (fConn == null) setupDb();
	}

	/* (non-Javadoc)
	 * @see java.util.Iterator#hasNext()
	 */
	public boolean hasNext() 
	{
		boolean out = false;
		try
		{
			if (curRS == null)
			{
			}
			else
			{
				out = curRS.next();
				curRS.previous(); // Roll back
			}
		}
		catch (Exception e)
		{
			Logger.getLogger(Logger.GLOBAL_LOGGER_NAME).severe("hasNext() error: " + e.getMessage());
			System.exit(1);
		}
		
		//logger.log(Level.INFO, "hasNext(): " + out);
		
		return out;
	}


	/* (non-Javadoc)
	 * @see java.util.Iterator#next()
	 */
	public Cpg next() {

		Cpg out = null;
		
		try
		{
			if (this.params.useReadTableSchema())
			{
				int onRead = 0;
				boolean finished = false;
				int curChrPos = -1;
				while (!finished)
				{
					finished = (!curRS.next());
					if (!finished)
					{
						int thisChromPos = curRS.getInt("cpg.chromPos");

						if (onRead == 0)
						{
							curChrPos = thisChromPos;
							out = new Cpg(thisChromPos,
									curRS.getString("cpg.strand").equals("-"),
									(short)0,
									(short)0,
									(short)0,
									(short)0,
									(short)0,
									(short)curRS.getInt("cpg.totalReadsOpposite"),
									(short)curRS.getInt("cpg.aReadsOpposite"),
									curRS.getInt("cpg.cpgWeight"),
									(short)0,
									(short)0,
									curRS.getString("cpg.nextBaseUpperCase").charAt(0)
							);
							//System.err.println("New CpG at " + thisChromPos + " strand=" + out.getStrandStr() + "\tdbstrand=" + curRS.getString("cpg.strand"));
						}

						// Are we finished?
						if (curChrPos != thisChromPos)
						{
							// System.err.println("\tFinished CpG at " + thisChromPos);
							finished = true;
							// Rewind
							curRS.previous(); // Roll back
						}
						else
						{
							//System.err.println("\tAdding to CpG at " + thisChromPos);
							CpgRead read = new CpgRead(
									curRS.getInt("cpg.readId"),
									(short)curRS.getInt("cpg.cRead"),
									(short)curRS.getInt("cpg.cReadNonconversionFilt"),
									(short)curRS.getInt("cpg.tRead"),
									(short)curRS.getInt("cpg.agRead"),
									(short)curRS.getInt("cpg.nextBaseGread"),
									curRS.getString("cpg.nextBaseUpperCase").charAt(0)
							);
							out.addRead(read);
							// Increment
							onRead++;
						}
					}
				}

				
			}
			else
			{
				curRS.next();

				out = new Cpg(curRS.getInt("cpg.chromPos"),
						curRS.getString("cpg.strand").equals("-"),
						(short)curRS.getInt("cpg.totalReads"),
						(short)curRS.getInt("cpg.cReads"),
						(short)curRS.getInt("cpg.cReadsNonconversionFilt"),
						(short)curRS.getInt("cpg.tReads"),
						(short)curRS.getInt("cpg.agReads"),
						(short)curRS.getInt("cpg.totalReadsOpposite"),
						(short)curRS.getInt("cpg.aReadsOpposite"),
						curRS.getInt("cpg.cpgWeight"),
						(short)curRS.getInt("cpg.nextBaseGreads"),
						(short)curRS.getInt("cpg.nextBaseTotalReads"),
						'0'
						);
			}
			
			// Downsample data.
			double downsamplingFactor = this.params.getSamplingFactor();
			if (downsamplingFactor != 0.0)
			{
				out = out.downsample(downsamplingFactor);
			}

			
		}
		catch (Exception e)
		{
			Logger.getLogger(Logger.GLOBAL_LOGGER_NAME).severe("next() error: " + e.getMessage());
			System.exit(1);
		}
		
		return out;
	}

	/* (non-Javadoc)
	 * @see java.util.Iterator#remove()
	 */
	public void remove() {
	}
	
	protected static String getSql(MethylDbQuerier params)
	throws Exception
	{
		return sqlHelper(params, null);
	}
	
	
	private static void fillPrep(MethylDbQuerier params, PreparedStatement prep)
	throws Exception
	{
		sqlHelper(params, prep);
	}
	

	protected static PreparedStatement getPrep(MethylDbQuerier inParams)
	throws Exception
	{
		String sql = getSql(inParams);
		return getPrep(sql);
	}
	
	protected static PreparedStatement getPrep(String sql)
	throws SQLException
	{
		PreparedStatement prep = cPreps.get(sql);
		if (prep==null)
		{
			prep = cConn.prepareStatement(sql);
			cPreps.put(sql, prep);
			Logger.getLogger(Logger.GLOBAL_LOGGER_NAME).severe("Making prepared statement: " + sql );
		}
		return prep;
	}

	
	/**
	 * @param params
	 * @param prep If not null, we actually fill in the prepared statment with values
	 * @return
	 */
	private static String sqlHelper(MethylDbQuerier params, PreparedStatement prep)
	throws Exception
	{
		// If there is a feature join, usually the feature table is MUCH smaller, so 
		// we put it first and do a straight_join so that it uses the feature index
		// first.  Here is an example
//		mysql> select count(*) from features_chr11, methylCGsRich_normal123009_chr11 cpg  WHERE ((cpg.chromPos>=0 AND cpg.chromPos<=5531960706)) AND  ((cpg.chromPos>=(chromPosStart)) AND (cpg.chromPos<=(chromPosEnd))  AND (featType = 'hg18.ES.H3K27me3.HMM.gtf'))  AND ((cpg.cReads+cpg.tReads) >= 4) AND ((cpg.totalReadsOpposite=0) OR ((cpg.aReadsOpposite/cpg.totalReadsOpposite)<=0.2)) AND (chromPosEnd>0) AND (chromPosStart<5531960706) ORDER BY chromPos ;
//		+----------+
//		| count(*) |
//		+----------+
//		|    49327 |
//		+----------+
//		1 row in set (6 min 7.95 sec)
//
//		mysql> select straight_join count(*) from features_chr11, methylCGsRich_normal123009_chr11 cpg  WHERE ((cpg.chromPos>=0 AND cpg.chromPos<=5531960706)) AND  ((cpg.chromPos>=(chromPosStart)) AND (cpg.chromPos<=(chromPosEnd))  AND (featType = 'hg18.ES.H3K27me3.HMM.gtf'))  AND ((cpg.cReads+cpg.tReads) >= 4) AND ((cpg.totalReadsOpposite=0) OR ((cpg.aReadsOpposite/cpg.totalReadsOpposite)<=0.2)) AND (chromPosEnd>0) AND (chromPosStart<5531960706) ORDER BY chromPos ;
//		+----------+
//		| count(*) |
//		+----------+
//		|    49327 |
//		+----------+
//		1 row in set (1 min 51.58 sec)

		
		String featTabSec = "";
		for (int i = 0; i < params.getFeatFilterTypes().size(); i++)
		{
			featTabSec += params.getFeatTable() + " feat" + i;
			featTabSec += ", ";
		}
		String joinSec = (params.usesFeatTable()) ? "straight_join" : "";
		
		String methTable = params.getMethylTable();
		String sql = String.format("select %s * from %s %s cpg WHERE ", joinSec, featTabSec, methTable);
		MethylDbQuerier.HelperOutput output = params.sqlWhereSecHelper(prep, "cpg");
		sql += output.sql;
		sql += (params.useReadTableSchema()) ? " GROUP BY chromPos,readId " : " GROUP BY chromPos "; // If you don't do this, you get multiple instances of the same CpG if it overlaps multiple features.
		sql += " ORDER BY chromPos ;";
		return sql;
	}
	
	
	protected static void setupDb()
	throws Exception
	{
		if (CpgIterator.cConn == null)
		{
			String connStr = MethylDbQuerier.connStr;
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

	
	
}


