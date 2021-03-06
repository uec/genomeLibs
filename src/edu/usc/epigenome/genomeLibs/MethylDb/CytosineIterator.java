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
 * @author yaping liu
 *
 */
public class CytosineIterator implements Iterator<Cytosine> {

	// Class vars
	protected static Connection cConn = null; 
	protected static Map<String,PreparedStatement> cPreps = new HashMap<String,PreparedStatement>();
	
	// Object vars
	protected MethylDbQuerier params = null;
	protected ResultSet curRS = null;
	protected int curNumRows = -1;
	protected String curColStrings = null;
	//protected double methyDens = Double.NaN;
	

	/**
	 * 
	 */
	public CytosineIterator(MethylDbQuerier inParams)
	throws Exception 
	{
		super();
		this.init(inParams);
	}

	public CytosineIterator(MethylDbQuerier inParams, String connStrCytosine, String sqlStatement)
	throws Exception 
	{
		super();
		this.init(inParams, connStrCytosine, sqlStatement);
	}
	
	//public double methyDensIterater(){
		//CytosineIterator it = new  
	//}
	
	//public double getMethyDens(){
	//	return methyDens;
	//}
	
	public int getCurNumRows() {
		return curNumRows;
	}
	
	public String getString (String columnName)
	throws Exception
	{
		String colStrings = curRS.getString(columnName);
		this.curColStrings = colStrings;
		return colStrings;
	}




	public CytosineIterator(String chrom, int startCoord, int endCoord, String inTablePrefix)
	throws Exception 
	{
		super();
		MethylDbQuerier params = new MethylDbQuerier();
		if (inTablePrefix!=null) params.methylTablePrefix = inTablePrefix;
		params.addRangeFilter(chrom, startCoord, endCoord);
		this.init(params);
	}
		
	
	public CytosineIterator(String chrom, String inTablePrefix)
	throws Exception 
	{
		super();
		MethylDbQuerier params = new MethylDbQuerier();
		if (inTablePrefix!=null) params.methylTablePrefix = inTablePrefix;
		params.addRangeFilter(chrom);
		this.init(params);

	}

	public int init(MethylDbQuerier inParams)
	throws Exception
	{
		this.params = inParams;

		// Check if we've started DB connection
		if (cConn == null) setupDb();
		
		PreparedStatement prep = CytosineIterator.getPrep(inParams);
		CytosineIterator.fillPrep(params, prep);
		Logger.getLogger(Logger.GLOBAL_LOGGER_NAME).fine("Starting query execute");
		curRS = prep.executeQuery();		
		curRS.last();
		int numRows = curRS.getRow();
		curRS.beforeFirst();
		Logger.getLogger(Logger.GLOBAL_LOGGER_NAME).fine("Finished query execute");
		//this.methyDens = methyDensIterater();
		this.curNumRows = numRows;
		return numRows;
	}
	
	public int init(MethylDbQuerier inParams, String connStrCytosine, String sqlStatement)
	throws Exception
	{
		this.params = inParams;

		// Check if we've started DB connection
		if (cConn == null) setupDb(connStrCytosine);
		
		PreparedStatement prep = CytosineIterator.getPrep(sqlStatement);
		//CytosineIterator.fillPrep(params, prep, true);
		Logger.getLogger(Logger.GLOBAL_LOGGER_NAME).fine("Starting query execute");
		curRS = prep.executeQuery();		
		curRS.last();
		int numRows = curRS.getRow();
		curRS.beforeFirst();
		Logger.getLogger(Logger.GLOBAL_LOGGER_NAME).fine("Finished query execute");
		//this.methyDens = methyDensIterater();
		this.curNumRows = numRows;
		return numRows;
	}
	
	public void init(MethylDbQuerier inParams, String connStrCytosine, String sqlStatement, boolean insertion)
	throws Exception
	{
		this.params = inParams;

		// Check if we've started DB connection
		if (cConn == null) setupDb(connStrCytosine);
		
		PreparedStatement prep = CytosineIterator.getPrep(sqlStatement);
		//CytosineIterator.fillPrep(params, prep, true);
		Logger.getLogger(Logger.GLOBAL_LOGGER_NAME).fine("Starting insertion execute");
		prep.executeQuery();		
		Logger.getLogger(Logger.GLOBAL_LOGGER_NAME).fine("Finished insertion execute");
		//this.methyDens = methyDensIterater();

	}
	
	public CytosineIterator()
	throws Exception 
	{
		throw new Exception("Whole genome Cytosine iterator not yet supported");
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
	public Cytosine next() {

		Cytosine out = null;
		try
		{
			curRS.next();
			
			out = new Cytosine(curRS.getInt("chromPos"),
					curRS.getString("strand").equals("-"),
					curRS.getShort("totalReads"),
					curRS.getShort("cReads"),
					curRS.getShort("cReadsNonconversionFilt"),
					curRS.getShort("tReads"),
					curRS.getShort("agReads"),
					curRS.getShort("totalReadsOpposite"),
					curRS.getShort("aReadsOpposite"),
					curRS.getString("preBaseRefUpperCase"),
					curRS.getString("nextBaseRefUpperCase"),
					curRS.getDouble("fracMeth"),
					curRS.getInt("gchWeight"),
					curRS.getInt("gcgWeight"),
					curRS.getInt("hcgWeight"));
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
	
	public Cytosine next(boolean asmFlag) {

		Cytosine out = null;
		try
		{
			curRS.next();
			
			out = new Cytosine(curRS.getInt("chromPos"),
					curRS.getString("strand").equals("-"),
					curRS.getShort("totalReads"),
					curRS.getShort("cReads"),
					curRS.getShort("cReadsNonconversionFilt"),
					curRS.getShort("tReads"),
					curRS.getShort("agReads"),
					curRS.getShort("totalReadsOpposite"),
					curRS.getShort("aReadsOpposite"),
					curRS.getInt("alleleChromPos"),
					curRS.getString("ABaseRefUpperCase"),
					curRS.getString("BBaseRefUpperCase"),
					curRS.getShort("ACReads"),
					curRS.getShort("BCReads"),
					curRS.getShort("ATReads"),
					curRS.getShort("BTReads"),
					curRS.getString("preBaseRefUpperCase"),
					curRS.getString("nextBaseRefUpperCase"));
					//curRS.getDouble("fracMeth"),
					//curRS.getDouble("fracAMeth"),
					//curRS.getDouble("fracBMeth"),
					//curRS.getInt("gchWeight"),
					//curRS.getInt("gcgWeight"),
					//curRS.getInt("hcgWeight"));
			// Downsample data.
			/*double downsamplingFactor = this.params.getSamplingFactor();
			if (downsamplingFactor != 0.0)
			{
				out = out.downsample(downsamplingFactor);
			}*/

			
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

		
		//String featTabSec = "";
		//for (int i = 0; i < params.getFeatFilterTypes().size(); i++)
		//{
			//featTabSec += params.getFeatTable() + " feat" + i;
			//featTabSec += ", ";
		//}
		//String joinSec = (params.usesFeatTable()) ? "straight_join" : "";
		
		String methTable = params.getMethylTable();
		//String methTable = params.methylTablePrefix;
		String sql = String.format("select * from %s WHERE ", methTable);
		MethylDbQuerier.HelperOutput output = params.sqlWhereSecHelperForGpc(prep, null);
		
		sql += output.sql;
		//sql += " GROUP BY chromPos "; // If you don't do this, you get multiple instances of the same CpG if it overlaps multiple features.
		sql += " ORDER BY chromPos ;";
		
		return sql;
	}
	
	
/*	private static String sqlHelperMethyFrac(MethylDbQuerier params, PreparedStatement prep)
	throws Exception
	{
		
		String methTable = params.getMethylTable();
		//String methTable = params.methylTablePrefix;
		String sql = String.format("select AVG(`fracMeth`) from %s WHERE ", methTable);
		MethylDbQuerier.HelperOutput output = params.sqlWhereSecHelperForGpc(prep, null);
		
		sql += output.sql;
		//sql += " GROUP BY chromPos "; // If you don't do this, you get multiple instances of the same CpG if it overlaps multiple features.
		sql += " ORDER BY chromPos ;";
		
		return sql;
	}
	*/
	protected static void setupDb()
	throws Exception
	{
		if (CytosineIterator.cConn == null)
		{
			String connStr = "jdbc:mysql://localhost/Cytosine";
			//String connStr = MethylDbQuerier.connStr;
			Class.forName("com.mysql.jdbc.Driver").newInstance();
			System.err.println("Getting connection for " + connStr);
			cConn = DriverManager.getConnection(connStr, "yaping", "lyping1986");
		}
		
	}
	
	protected static void setupDb(String connStrCytosine)
	throws Exception
	{
		if (CytosineIterator.cConn == null)
		{
			String connStr = connStrCytosine;
			//String connStr = MethylDbQuerier.connStr;
			Class.forName("com.mysql.jdbc.Driver").newInstance();
			System.err.println("Getting connection for " + connStr);
			cConn = DriverManager.getConnection(connStrCytosine, "yaping", "lyping1986");
		}
		
	}
	
	protected static void cleanupDb()
	throws Exception
	{
		cConn.close();
	}
	
	// boolean absolute( int row ) throws SQLException;
	
	
}
