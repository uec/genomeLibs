package edu.usc.epigenome.genomeLibs.FeatDb;

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

import org.biojava.bio.program.gff.GFFRecord;
import org.biojava.bio.program.gff.SimpleGFFRecord;
import org.biojava.bio.seq.StrandedFeature;
import org.biojava.bio.seq.StrandedFeature.Strand;
import org.usckeck.genome.GFFUtils;


/**
 * @author benb
 *
 */
public class FeatIterator implements Iterator<GFFRecord> {

	// Class vars
	protected static Connection cConn = null; 
	protected static Map<String,PreparedStatement> cPreps = new HashMap<String,PreparedStatement>();
	private static Logger logger = Logger.getLogger(Logger.GLOBAL_LOGGER_NAME); // "edu.usc.epigenome.genomeLibs.MethylDb.CpgIterator");
	
	// Object vars
	protected FeatDbParams params = null;
	protected ResultSet curRS = null;
	

	/**
	 * 
	 */
	public FeatIterator(FeatDbParams inParams)
	throws Exception 
	{
		super();
		this.init(inParams);
	}



	
	public FeatIterator(String chrom, String inTablePrefix, int startCoord, int endCoord, String featType)
	throws Exception 
	{
		super();
		FeatDbParams params = new FeatDbParams();
		if (inTablePrefix!=null) params.tablePrefix = inTablePrefix;
		params.addRangeFilter(chrom, startCoord, endCoord);
		params.addFeatFilter(featType);
		this.init(params);
	}
	
	public FeatIterator(String chrom, String inTablePrefix, int startCoord, int endCoord)
	throws Exception 
	{
		super();
		FeatDbParams params = new FeatDbParams();
		if (inTablePrefix!=null) params.tablePrefix = inTablePrefix;
		params.addRangeFilter(chrom, startCoord, endCoord);
		this.init(params);
	}
		
	
	public FeatIterator(String chrom, String inTablePrefix)
	throws Exception 
	{
		super();
		FeatDbParams params = new FeatDbParams();
		if (inTablePrefix!=null) params.tablePrefix = inTablePrefix;
		params.addRangeFilter(chrom);
		this.init(params);

	}

	public void init(FeatDbParams inParams)
	throws Exception
	{
		this.params = inParams;

		// Check if we've started DB connection
		if (cConn == null) setupDb();
		PreparedStatement prep = FeatIterator.getPrep(inParams);
		FeatIterator.fillPrep(params, prep);
		curRS = prep.executeQuery();		
	}
	
	
	
	public FeatIterator()
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
	public GFFRecord next() {

		GFFRecord out = null;
		try
		{
			curRS.next();
			
			String strandStr = curRS.getString("strand");
			Strand strand = StrandedFeature.UNKNOWN;
			if (strandStr.equals("+"))
			{
				strand = StrandedFeature.POSITIVE;
			}
			else if (strandStr.equals("-"))
			{
				strand = StrandedFeature.NEGATIVE;
			}
				
			
			out = new SimpleGFFRecord(
			"chr",
            curRS.getString("featType"),
			"exon",
            curRS.getInt("chromPosStart"),
            curRS.getInt("chromPosEnd"),
            curRS.getDouble("score"),
            strand,
            0,
            "",
            new HashMap(1));
	
			String name = curRS.getString("name");
			if (name != null)
			{
				GFFUtils.setGffRecordName(out, name);
			}
		}
		catch (Exception e)
		{
			logger.log(Level.SEVERE, "next() error: " + e.toString());
			e.printStackTrace();
			System.exit(1);
		}
		
		return out;
	}

	/* (non-Javadoc)
	 * @see java.util.Iterator#remove()
	 */
	public void remove() {
	}
	
	protected static String getSql(FeatDbParams params)
	throws Exception
	{
		return sqlHelper(params, null);
	}
	
	
	private static void fillPrep(FeatDbParams params, PreparedStatement prep)
	throws Exception
	{
		sqlHelper(params, prep);
	}
	

	protected static PreparedStatement getPrep(FeatDbParams inParams)
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
	private static String sqlHelper(FeatDbParams params, PreparedStatement prep)
	throws Exception
	{
		String table = params.getTable();
		String sql = String.format("select * from %s feat WHERE ", table);
		sql += params.sqlWhereSecHelper(prep, "feat");
		sql += " ORDER BY chromPosStart ;";
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
