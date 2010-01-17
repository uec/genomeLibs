package edu.usc.epigenome.genomeLibs.MethylDb;

import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.PreparedStatement;
import java.sql.ResultSet;

import java.sql.SQLException;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.logging.Level;
import java.util.logging.Logger;

import edu.usc.epigenome.genomeLibs.ListUtils;


/**
 * @author benb
 *
 */
public class CpgIteratorMultisample implements Iterator<Cpg[]> {

	// Class vars
	protected static Connection cConn = null; 
	protected static Map<String,PreparedStatement> cPreps = new HashMap<String,PreparedStatement>();
	
	// Object vars
	protected MethylDbQuerier params = null;
	protected List<String> sampleTablePrefixes = null;
	protected ResultSet curRS = null;
	protected int curNumRows = -1;
	protected boolean allowNumRows = true; // If we set this to false, we use A LOT less memory
	

	/**
	 * 
	 */
	public CpgIteratorMultisample(MethylDbQuerier inParams, List<String> inSampleTablePrefixes)
	throws Exception 
	{
		super();
		this.init(inParams, inSampleTablePrefixes);
	}

	/**
	 * @param inParams
	 * @param inSampleTablePrefixes
	 * @param inAllowNumRows If we set this to false, we use A LOT less memory, but can't use "getNumRows" in advance
	 * @throws Exception
	 */
	public CpgIteratorMultisample(MethylDbQuerier inParams, List<String> inSampleTablePrefixes, boolean inAllowNumRows)
	throws Exception 
	{
		super();
		this.allowNumRows = inAllowNumRows;
		this.init(inParams, inSampleTablePrefixes);
	}
	
	public int getCurNumRows() {
		return curNumRows;
	}



	public int init(MethylDbQuerier inParams, List<String> inSampleTablePrefixes)
	throws Exception
	{
		this.params = inParams;
		this.sampleTablePrefixes = inSampleTablePrefixes;
		
		// Check if we've started DB connection
		if (cConn == null) setupDb();
		PreparedStatement prep = this.getPrep(inParams);
		this.fillPrep(params, prep);
		Logger.getLogger(Logger.GLOBAL_LOGGER_NAME).fine("Starting query execute");
//		if (!this.allowNumRows)
//		{
//			System.err.println("Setting Featch size to 8192");
//			prep.setFetchSize(8192);
//		}
		curRS = prep.executeQuery();		
		Logger.getLogger(Logger.GLOBAL_LOGGER_NAME).fine("Finished query execute");
		int numRows = -1;
		if (this.allowNumRows)
		{
			Logger.getLogger(Logger.GLOBAL_LOGGER_NAME).fine("About to advance to end to count Cpgs");
			curRS.last();
			numRows = curRS.getRow();
			curRS.beforeFirst();
			this.curNumRows = numRows;
		}
		Logger.getLogger(Logger.GLOBAL_LOGGER_NAME).fine("Finished query");
		return numRows;
	}
	
	
	
	public CpgIteratorMultisample()
	throws Exception 
	{
		throw new Exception("Whole genome Cpg iterator not yet supported");
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
	public Cpg[] next() {

		int nS = this.sampleTablePrefixes.size();
		Cpg[] out = new Cpg[nS];
		try
		{
			curRS.next();
			
			for (int i = 0; i < nS; i++)
			{
				Cpg cpg = new Cpg(curRS.getInt("cpg" + i + ".chromPos"),
						curRS.getString("cpg" + i + ".strand").equals("-"),
						curRS.getShort("cpg" + i + ".totalReads"),
						curRS.getShort("cpg" + i + ".cReads"),
						curRS.getShort("cpg" + i + ".cReadsNonconversionFilt"),
						curRS.getShort("cpg" + i + ".tReads"),
						curRS.getShort("cpg" + i + ".agReads"),
						curRS.getShort("cpg" + i + ".totalReadsOpposite"),
						curRS.getShort("cpg" + i + ".aReadsOpposite"));

				// Downsample data.
				double downsamplingFactor = this.params.getSamplingFactor();
				if (downsamplingFactor != 0.0)
				{
					cpg = cpg.downsample(downsamplingFactor);
				}

				out[i] = cpg;
				//System.err.printf("out[%d] = %s\n",i,cpg.toString());
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
	
	protected String getSql(MethylDbQuerier params)
	throws Exception
	{
		return sqlHelper(params, null);
	}
	
	
	private void fillPrep(MethylDbQuerier params, PreparedStatement prep)
	throws Exception
	{
		sqlHelper(params, prep);
	}
	

	protected  PreparedStatement getPrep(MethylDbQuerier inParams)
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
	private String sqlHelper(MethylDbQuerier params, PreparedStatement prep)
	throws Exception
	{
		int nS = this.sampleTablePrefixes.size();

		
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

		
		String featTabSec = (params.usesFeatTable()) ? (params.getFeatTable() + ", " ) : "";
		String joinSec = (params.usesFeatTable()) ? " straight_join " : "";
		
		// Table sec
		String chrom = params.getChrom();
		String[] selectSecs = new String[nS];
		String[] tableSecs = new String[nS];
		for (int i = 0; i < nS; i++)
		{
			selectSecs[i] = String.format("cpg%d.*", i);
			tableSecs[i] = String.format("%s cpg%d", this.sampleTablePrefixes.get(i) + chrom, i);
		}
		ListUtils.setDelim(", ");
		String sql = "SELECT ";
		sql += joinSec;
		sql += ListUtils.excelLine(selectSecs);
		sql += " FROM ";
		sql += featTabSec;
		sql += ListUtils.excelLine(tableSecs);
		sql += " WHERE ";

		// Constraints are true for all samples
		String[] secs = new String[nS];
		int curInd = 1;
		for (int i = 0; i < nS; i++)
		{
			MethylDbQuerier.HelperOutput output = params.sqlWhereSecHelper(prep, "cpg"+i, curInd); 
			secs[i] = output.sql;
			int newCurInd = output.newCurInd;
			
			// Weird way to get the number of params.  Stupid Regex functions in eclipse at
			// least can't use "/?" sequences
			int lastFound = -1;
			int numParams = 0;
			while ((lastFound = secs[i].indexOf('?', lastFound+1)) > 0)
			{
				numParams++;
			}
			curInd += numParams;
			if ((prep!= null) && (newCurInd != curInd)) throw new Exception ("Found " + curInd + " params with old method, " + newCurInd + " with new method: "
					+ output.sql);
			curInd = newCurInd;
		}		
		ListUtils.setDelim(" AND ");
		sql += ListUtils.excelLine(secs);
		
		
		// Join the samples
		secs = new String[nS-1];
		for (int i = 1; i < nS; i++)
		{
			secs[i-1] = String.format("(cpg0.chromPos=cpg%d.chromPos)",i);
		}		
		ListUtils.setDelim(" AND ");
		if (secs.length > 0) sql += " AND ";
		sql += ListUtils.excelLine(secs);
		
		// And finish
		sql += " ORDER BY cpg0.chromPos;";
		return sql;
	}
	
	
	protected void setupDb()
	throws Exception
	{
		String connStr = this.params.connStr;
		Class.forName("com.mysql.jdbc.Driver").newInstance();
		System.err.println("Getting connection for " + connStr);
		cConn = DriverManager.getConnection(connStr);
		
	}
	
	protected static void cleanupDb()
	throws Exception
	{
		cConn.close();
	}

	
	
}
