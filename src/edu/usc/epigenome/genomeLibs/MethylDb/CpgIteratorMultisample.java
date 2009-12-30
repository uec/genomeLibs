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
	private static Logger logger = Logger.getLogger(Logger.GLOBAL_LOGGER_NAME); // "edu.usc.epigenome.genomeLibs.MethylDb.CpgIterator");
	
	// Object vars
	protected MethylDbParams params = null;
	protected List<String> sampleTablePrefixes = null;
	protected ResultSet curRS = null;
	

	/**
	 * 
	 */
	public CpgIteratorMultisample(MethylDbParams inParams, List<String> inSampleTablePrefixes)
	throws Exception 
	{
		super();
		this.init(inParams, inSampleTablePrefixes);
	}




	public void init(MethylDbParams inParams, List<String> inSampleTablePrefixes)
	throws Exception
	{
		this.params = inParams;
		this.sampleTablePrefixes = inSampleTablePrefixes;

		// Check if we've started DB connection
		if (cConn == null) setupDb();
		PreparedStatement prep = this.getPrep(inParams);
		this.fillPrep(params, prep);
		curRS = prep.executeQuery();		
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
			logger.log(Level.SEVERE, "hasNext() error: " + e.getMessage());
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
				
				out[i] = cpg;
				//System.err.printf("out[%d] = %s\n",i,cpg.toString());
			}
			
		}
		catch (Exception e)
		{
			logger.log(Level.SEVERE, "next() error: " + e.getMessage());
			System.exit(1);
		}
		
		return out;
	}

	/* (non-Javadoc)
	 * @see java.util.Iterator#remove()
	 */
	public void remove() {
	}
	
	protected String getSql(MethylDbParams params)
	throws Exception
	{
		return sqlHelper(params, null);
	}
	
	
	private void fillPrep(MethylDbParams params, PreparedStatement prep)
	throws Exception
	{
		sqlHelper(params, prep);
	}
	

	protected  PreparedStatement getPrep(MethylDbParams inParams)
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
			logger.log(Level.SEVERE, "Making prepared statement: " + sql );
		}
		return prep;
	}

	
	/**
	 * @param params
	 * @param prep If not null, we actually fill in the prepared statment with values
	 * @return
	 */
	private String sqlHelper(MethylDbParams params, PreparedStatement prep)
	throws Exception
	{
		int nS = this.sampleTablePrefixes.size();

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
		sql += ListUtils.excelLine(selectSecs);
		sql += " FROM ";
		sql += ListUtils.excelLine(tableSecs);
		sql += " WHERE ";

		// Constraints are true for all samples
		String[] secs = new String[nS];
		int curInd = 1;
		for (int i = 0; i < nS; i++)
		{
			secs[i] = params.sqlWhereSecHelper(prep, "cpg"+i, curInd);
			
			// Weird way to get the number of params.  Stupid Regex functions in eclipse at
			// least can't use "/?" sequences
			int lastFound = -1;
			int numParams = 0;
			while ((lastFound = secs[i].indexOf('?', lastFound+1)) > 0)
			{
				numParams++;
			}
			//System.err.println("Found " + numParams + " params");
			curInd += numParams;
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
		sql += " ORDER BY cpg1.chromPos;";
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