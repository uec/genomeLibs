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
	private static Logger logger = Logger.getLogger(Logger.GLOBAL_LOGGER_NAME); // "edu.usc.epigenome.genomeLibs.MethylDb.CpgIterator");
	
	// Object vars
	protected MethylDbParams params = null;
	protected ResultSet curRS = null;
	

	/**
	 * 
	 */
	public CpgIterator(MethylDbParams inParams)
	throws Exception 
	{
		super();
		this.init(inParams);
	}



	
	public CpgIterator(String chrom, int startCoord, int endCoord, String inTablePrefix)
	throws Exception 
	{
		super();
		MethylDbParams params = new MethylDbParams();
		if (inTablePrefix!=null) params.tablePrefix = inTablePrefix;
		params.addRangeFilter(chrom, startCoord, endCoord);
		this.init(params);
	}
		
	
	public CpgIterator(String chrom, String inTablePrefix)
	throws Exception 
	{
		super();
		MethylDbParams params = new MethylDbParams();
		if (inTablePrefix!=null) params.tablePrefix = inTablePrefix;
		params.addRangeFilter(chrom);
		this.init(params);

	}

	public void init(MethylDbParams inParams)
	throws Exception
	{
		this.params = inParams;

		// Check if we've started DB connection
		if (cConn == null) setupDb();
		PreparedStatement prep = CpgIterator.getPrep(inParams);
		CpgIterator.fillPrep(params, prep);
		curRS = prep.executeQuery();		
	}
	
	
	
	public CpgIterator()
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
	public Cpg next() {

		Cpg out = null;
		try
		{
			curRS.next();
			
			out = new Cpg(curRS.getInt("chromPos"),
					curRS.getString("strand").equals("-"),
					curRS.getShort("totalReads"),
					curRS.getShort("cReads"),
					curRS.getShort("cReadsNonconversionFilt"),
					curRS.getShort("tReads"),
					curRS.getShort("agReads"),
					curRS.getShort("totalReadsOpposite"),
					curRS.getShort("aReadsOpposite"));
			
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
	
	protected static String getSql(MethylDbParams params)
	throws Exception
	{
		return sqlHelper(params, null);
	}
	
	
	private static void fillPrep(MethylDbParams params, PreparedStatement prep)
	throws Exception
	{
		sqlHelper(params, prep);
	}
	

	protected static PreparedStatement getPrep(MethylDbParams inParams)
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
	private static String sqlHelper(MethylDbParams params, PreparedStatement prep)
	throws Exception
	{
		String table = params.getTable();
		String sql = String.format("select * from %s cpg WHERE ", table);
		sql += params.sqlWhereSecHelper(prep, "cpg");
		sql += " ORDER BY chromPos ;";
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
