package edu.usc.epigenome.genomeLibs.FeatDb;

import java.sql.PreparedStatement;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.logging.Level;
import java.util.logging.Logger;

import edu.usc.epigenome.genomeLibs.ListUtils;
import edu.usc.epigenome.genomeLibs.GenomicRange.GenomicRange;

public class FeatDbQuerier {

	public final static String DEFAULT_TABLE_PREFIX = "features_";
	public final static String DEFAULT_CONN_STR = "jdbc:mysql://localhost/cr?user=benb";
	
	public final static String connStr = DEFAULT_CONN_STR;
	
	protected String tablePrefix = DEFAULT_TABLE_PREFIX;
	
	
	// Ranges
	protected Set<GenomicRange> rangeFilters = new HashSet<GenomicRange>(1);
	protected Set<String> featFilters = new HashSet<String>(10);

	
	public void addFeatFilter(String featType)
	{
		featFilters.add(featType);
	}
	
	public void addRangeFilter(String chrom)
	{
		GenomicRange gr = GenomicRange.fullChromRange(chrom, "hg18");
		this.addRangeFilter(gr);
	}
	
	public void addRangeFilter(String chrom, int startCoord, int endCoord)
	{
		GenomicRange gr = new GenomicRange(chrom, startCoord, endCoord);
		this.addRangeFilter(gr);
	}
	
	public void addRangeFilter(GenomicRange inGr)
	{
		this.rangeFilters.add(inGr);
	}
	
	public Set<String> getRangeFilterChroms()
	{
		Set<String> out = new HashSet<String>(rangeFilters.size());
		for (GenomicRange range : this.rangeFilters)
		{
			out.add(range.getChrom());
		}
		return out;
	}

	public String getChrom()
	throws Exception
	{
		Set<String> chroms = this.getRangeFilterChroms();
		String chrom;
		if (chroms.size() == 0)
		{
			Exception e = new Exception("can not get CpGs without a range filter");
			Logger.getLogger(Logger.GLOBAL_LOGGER_NAME).throwing("CpgIterator", "sqlHelper", e);
			throw e;
		}
		else if (chroms.size() > 1)
		{
			Exception e = new Exception("can not get CpGs from multiple chroms with a single SQL statement");
			Logger.getLogger(Logger.GLOBAL_LOGGER_NAME).throwing("CpgIterator", "sqlHelper", e);
			throw e;
		}
		else
		{
			chrom = chroms.toArray(new String[1])[0];
		}
		return chrom;
	}
	
	public String getTable()
	throws Exception
	{
		String chrom = this.getChrom();
		String table = this.tablePrefix + chrom;

		return table;
	}

	/********** GETTERS AND SETTERS ********/
	
	public String getTablePrefix() {
		return tablePrefix;
	}

	public void setTablePrefix(String tablePrefix) {
		this.tablePrefix = tablePrefix;
	}

	
	/******* PUBLIC DB STUFF ********/
	
	
	
	/**
	 * @param asName This is the table alias, i.e. SELECT ... FROM tab1 t1, tab2 t2 ... 
	 * @return
	 */
	public String sqlWhereSec(String asName)
	{
		String out = null;
		try
		{
			// I don't think this can really throw an exception if it's not 
			// setting a prepared statement.
			FeatDbQuerier.HelperOutput output = sqlWhereSecHelper(null, asName);
			out = output.sql;
		}
		catch (Exception e)
		{
			Logger.getLogger(Logger.GLOBAL_LOGGER_NAME).severe("Shouldn't have gotten here");
		}
		return out;
	}

	/**
	 * @param prep If supplied, we not only return the SQL but fill in the prep with its correct values
	 * @param asName the <alias name> of the table in the FROM section of the query
	 * @return the new current index
	 * @throws SQLException
	 */
	public int sqlWhereSecFillPrep(PreparedStatement prep, String asName)
	throws SQLException
	{
		FeatDbQuerier.HelperOutput output = sqlWhereSecHelper(prep, asName);
		return output.newCurInd;
	}

	/**
	 * @param prep If supplied, we not only return the SQL but fill in the prep with its correct values
	 * @param asName the <alias name> of the table in the FROM section of the query
	 * @return the new current index
	 * @throws SQLException
	 */
	public int sqlWhereSecFillPrep(PreparedStatement prep, String asName, int curInd)
	throws SQLException
	{
		FeatDbQuerier.HelperOutput output = sqlWhereSecHelper(prep, asName, curInd);
		return output.newCurInd;
	}

	
	// Private building of the where sec
	
	/**
	 * I'm making this public so that it can be used by MethylDbQuerier, but it generally should
	 * not be used.
	 * @param prep If supplied, we not only return the SQL but fill in the prep with its correct values
	 * @param asName the <alias name> of the table in the FROM section of the query
	 * @return
	 */
	public FeatDbQuerier.HelperOutput sqlWhereSecHelper(PreparedStatement prep, String asName)
	throws SQLException
	{
		return this.sqlWhereSecHelper(prep, asName, 1);
	}
	
	/**
	 * I'm making this public so that it can be used by MethylDbQuerier, but it generally should
	 * not be used.
	 * @param prep If supplied, we not only return the SQL but fill in the prep with its correct values
	 * @param asName the <alias name> of the table in the FROM section of the query
	 * @param curInd the offset of the current index number in the preparedStatement 
	 * @return
	 */
	public FeatDbQuerier.HelperOutput sqlWhereSecHelper(PreparedStatement prep, String asName, int curInd)
	throws SQLException
	{
		String asSec = (asName==null) ? "" : (asName+".");

		List<String> clauses = new ArrayList<String>(20);
		
		
		// Range filters
		if (this.rangeFilters.size()>0)
		{
			List<String> grClauses = new ArrayList<String>(this.rangeFilters.size());
			for (GenomicRange gr : this.rangeFilters)
			{
				// This is the best way to pick anything that OVERLAPS the range.
				// Notice than CONTAINED WITHIN the range would be a different query
				
				
				String clause = String.format("((?<=%schromPosEnd) AND (?>=%schromPosStart))",asSec, asSec);
				
				
				grClauses.add(clause);
				if (prep!=null)
				{
					prep.setInt(curInd++, gr.getStart());
					prep.setInt(curInd++, gr.getEnd());
				}

//   OLD WAY.  FOR SOME REASON I THOUGHT THIS WAS BETTER?
//				// This is the best way to pick anything that OVERLAPS the range.
//				// Notice than CONTAINED WITHIN the range would be a different query
//				String clause = String.format("(!(%schromPosStart<? AND %schromPosEnd<?) AND !(%schromPosStart>? AND %schromPosEnd>?))", asSec, asSec,asSec,asSec);
//				grClauses.add(clause);
//				if (prep!=null)
//				{
//					prep.setInt(curInd++, gr.getStart());
//					prep.setInt(curInd++, gr.getStart());
//					prep.setInt(curInd++, gr.getEnd());
//					prep.setInt(curInd++, gr.getEnd());
//				}
			}
			ListUtils.setDelim(" OR "); // Notice that it's a UNION of range filters
			clauses.add("(" + ListUtils.excelLine(grClauses.toArray(new String[1])) + ")");
		}
		
		// feat filters
		if (this.featFilters.size()>0)
		{
			List<String> fClauses = new ArrayList<String>(this.featFilters.size());
			for (String featFilter : this.featFilters)
			{
				String clause = String.format("(featType = ?)");
				fClauses.add(clause);
				if (prep!=null)
				{
					prep.setString(curInd++, featFilter);
				}
			}
			ListUtils.setDelim(" OR "); // Notice that it's a UNION of feat filters
			clauses.add("(" + ListUtils.excelLine(fClauses.toArray(new String[1])) + ")");
		}

		
		
		ListUtils.setDelim(" AND ");
		String sql = ListUtils.excelLine(clauses.toArray(new String[1]));
		
		return new FeatDbQuerier.HelperOutput(sql, curInd);
	}
	
	public class HelperOutput
	{
		/**
		 * @param sql
		 * @param newCurInd
		 */
		public HelperOutput(String sql, int newCurInd) {
			super();
			this.sql = sql;
			this.newCurInd = newCurInd;
		}
		public String sql = null;
		public int newCurInd = 0;
	}
	

	
	
}
