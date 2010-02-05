package edu.usc.epigenome.genomeLibs.FeatDb;

import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.Statement;

import java.sql.SQLException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.logging.Level;
import java.util.logging.Logger;

import org.biojava.bio.program.gff.GFFRecord;
import org.biojava.bio.program.gff.SimpleGFFRecord;
import org.biojava.bio.seq.StrandedFeature;
import org.biojava.bio.seq.StrandedFeature.Strand;
import org.usckeck.genome.GFFUtils;

import edu.usc.epigenome.genomeLibs.ListUtils;
import edu.usc.epigenome.genomeLibs.GenomicRange.GenomicRange;
import edu.usc.epigenome.genomeLibs.MethylDb.MethylDbUtils;


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
	protected FeatDbQuerier params = null;
	protected ResultSet curRS = null;
	

	/**
	 * 
	 */
	public FeatIterator(FeatDbQuerier inParams)
	throws Exception 
	{
		super();
		this.init(inParams);
	}



	
	public FeatIterator(String chrom, String inTablePrefix, int startCoord, int endCoord, String featType)
	throws Exception 
	{
		super();
		FeatDbQuerier params = new FeatDbQuerier();
		if (inTablePrefix!=null) params.tablePrefix = inTablePrefix;
		params.addRangeFilter(chrom, startCoord, endCoord);
		params.addFeatFilter(featType);
		this.init(params);
	}
	
	public FeatIterator(String chrom, String inTablePrefix, int startCoord, int endCoord)
	throws Exception 
	{
		super();
		FeatDbQuerier params = new FeatDbQuerier();
		if (inTablePrefix!=null) params.tablePrefix = inTablePrefix;
		params.addRangeFilter(chrom, startCoord, endCoord);
		this.init(params);
	}
		
	
	public FeatIterator(String chrom, String inTablePrefix)
	throws Exception 
	{
		super();
		FeatDbQuerier params = new FeatDbQuerier();
		if (inTablePrefix!=null) params.tablePrefix = inTablePrefix;
		params.addRangeFilter(chrom);
		this.init(params);

	}

	public void init(FeatDbQuerier inParams)
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
	
			String name = curRS.getString("refseqId");
			if ((name != null) && (!name.equalsIgnoreCase("null")))
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
	
	protected static String getSql(FeatDbQuerier params)
	throws Exception
	{
		return sqlHelper(params, null);
	}
	
	
	private static void fillPrep(FeatDbQuerier params, PreparedStatement prep)
	throws Exception
	{
		sqlHelper(params, prep);
	}
	

	protected static PreparedStatement getPrep(FeatDbQuerier inParams)
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
	private static String sqlHelper(FeatDbQuerier params, PreparedStatement prep)
	throws Exception
	{
		String table = params.getTable();
		String sql = String.format("select * from %s feat0 WHERE ", table);
		FeatDbQuerier.HelperOutput output = params.sqlWhereSecHelper(prep, "feat0");
		sql += output.sql;
		sql += " ORDER BY chromPosStart ;";
		return sql;
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

	/******
	 * STATIC UTILITY FUNCTIONS
	 */
	public static List<String> AllFeatTablePrefixes()
	throws Exception
	{
		return AllFeatTablePrefixes(FeatDbQuerier.DEFAULT_TABLE_PREFIX + "chr11");
	}
	
	public static List<String> AllFeatTablePrefixes(String table)
	throws Exception
	{
		setupDb();
		
		String sql = String.format("select distinct(featType) from %s;",table);

		Statement st = cConn.createStatement();
		ResultSet rs = st.executeQuery(sql);

		List<String> out = new ArrayList<String>(30);
		while (rs.next())
		{
			out.add(rs.getString(1));
		}
		
		rs.close();
		st.close();
		return out;
	}
	
	public static int TotalFeatSize(List<String> featTypes, Set<GenomicRange> rangeFilters, boolean intersection)
	throws Exception
	{
		int totalLen = 0;
		if ((rangeFilters==null) || (rangeFilters.size()==0))
		{
			// Special case
			totalLen = TotalFeatSize(featTypes, intersection);
		}
		else
		{
			FeatDbQuerier querier = new FeatDbQuerier();
			for (String featType : featTypes) { querier.addFeatFilter(featType); }
			for (GenomicRange range : rangeFilters) { querier.addRangeFilter(range); }
			FeatIterator it = new FeatIterator(querier);
			
			while (it.hasNext())
			{
				GFFRecord rec = it.next();
				totalLen += (rec.getEnd() - rec.getStart() + 1);
			}
		}
		
		return totalLen;
	}

	


	public static int TotalFeatSize(List<String> featTypes, boolean intersection)
	throws Exception
	{
		return TotalFeatSize(featTypes, intersection, MethylDbUtils.CHROMS);
	}
	
	public static int TotalFeatSize(List<String> featTypes, boolean intersection, List<String> chroms)
	throws Exception
	{
		int total = 0;
		for (String chr : chroms)
		{
			total += TotalFeatSize(featTypes, intersection, chr);
		}
		return total;
	}
	
	public static int TotalFeatSize(List<String> featTypes, boolean intersection, String chr)
	throws Exception
	{
		setupDb();
		
		String[] clauses = new String[featTypes.size()];
		for (int i = 0; i < featTypes.size(); i++)
		{
			clauses[i] = String.format("featType = '%s'", featTypes.get(i));
		}
		ListUtils.setDelim( (intersection) ? " AND " : " OR ");
		String whereClause = ListUtils.excelLine(clauses);
		
		String sql = String.format("SELECT SUM(chromPosEnd-chromPosStart) FROM %s WHERE %s;",
				FeatDbQuerier.DEFAULT_TABLE_PREFIX + chr, whereClause);
		Logger.getLogger(Logger.GLOBAL_LOGGER_NAME).severe(sql);

		Statement st = cConn.createStatement();
		ResultSet rs = st.executeQuery(sql);

		int out = 0;
		if (rs.next())
		{
			out = rs.getInt(1);
		}
		else
		{
			throw new Exception("Following SQL doesn't return a result:\n" + sql);
		}
		
		rs.close();
		st.close();
		return out;
	}


	
}
