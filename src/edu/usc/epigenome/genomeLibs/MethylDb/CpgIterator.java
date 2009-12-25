package edu.usc.epigenome.genomeLibs.MethylDb;

import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * @author benb
 *
 */
/**
 * @author benb
 *
 */
public class CpgIterator implements Iterator<Cpg> {

	final private static String TABLE_PREFIX = "methylCGsRich_tumor_";
//	final private static String TABLE_PREFIX = "RestingNucleosomes_CD4_";

	// Class vars
	protected static Connection cConn = null; 
	protected static Map<String,PreparedStatement> cByCoordsPreps = new HashMap<String,PreparedStatement>();
	protected static Map<String,PreparedStatement> cByChromPreps = new HashMap<String,PreparedStatement>();
	private static Logger logger = Logger.getLogger(Logger.GLOBAL_LOGGER_NAME); // "edu.usc.epigenome.genomeLibs.MethylDb.CpgIterator");
	
	// Object vars
	protected ResultSet curRS = null;
	

	/**
	 * 
	 */
	public CpgIterator(String chrom, int startCoord, int endCoord)
	throws Exception 
	{
		super();
		
		// Check if we've started DB connection
		if (cConn == null) setupDb();
		
		PreparedStatement prep = cByCoordsPreps.get(chrom);
		if (prep == null)
		{
			String table = TABLE_PREFIX + chrom;
			String sql = "select * from " + table + " WHERE chromPos >= ? AND chromPos <= ?;";
			prep = cConn.prepareStatement(sql);
			cByCoordsPreps.put(chrom, prep);
			logger.log(Level.INFO, "Making prepared statement for chrom " + chrom + ": " + sql );
		}
		else
		{
			//logger.log(Level.INFO, "Found prepared statement for chrom " + chrom);
		}

		prep.setInt(1, startCoord);
		prep.setInt(2, endCoord);
		curRS = prep.executeQuery();
		
	}
		
	public CpgIterator(String chrom)
	throws Exception 
	{
		super();
		
		// Check if we've started DB connection
		if (cConn == null) setupDb();
		
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
	
	
	
	protected static void setupDb()
	throws Exception
	{
		String connStr = "jdbc:mysql://localhost/cr?user=benb";
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
