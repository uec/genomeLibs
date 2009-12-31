package edu.usc.epigenome.genomeLibs.MethylDb;

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

public class MethylDbQuerier {

	public final static String DEFAULT_METHYL_TABLE_PREFIX = "methylCGsRich_tumor10x_";
	public final static String DEFAULT_FEAT_TABLE_PREFIX = "features_";
	public final static String DEFAULT_CONN_STR = "jdbc:mysql://localhost/cr?user=benb";
	
	
	public String methylTablePrefix = DEFAULT_METHYL_TABLE_PREFIX;
	public String featureTablePrefix = DEFAULT_FEAT_TABLE_PREFIX;
	public String connStr = DEFAULT_CONN_STR;
	
	
	// Ranges
	protected Set<GenomicRange> rangeFilters = new HashSet<GenomicRange>(1);
	protected Set<MethylDbQuerier.FeatClass> featFilters = new HashSet<MethylDbQuerier.FeatClass>(1);

	
	// Filtering
	protected int minCTreads = 0;
	protected boolean useNonconversionFilter = true;
	protected double maxOppstrandAfrac = 0.2; // Double.MAX_VALUE;
	
	
	
	public void addFeatFilter(String featType)
	{
		this.addFeatFilter(featType, 0);
	}
	
	
	public void addFeatFilter(String featType, int flankSize)
	{
		featFilters.add(new MethylDbQuerier.FeatClass(featType, flankSize));
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
		String table = this.methylTablePrefix + chrom;

		return table;
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
			out = sqlWhereSecHelper(null, asName);
		}
		catch (Exception e)
		{
			Logger.getLogger(Logger.GLOBAL_LOGGER_NAME).severe("Shouldn't have gotten here");
		}
		return out;
	}

	public void sqlWhereSecFillPrep(PreparedStatement prep, String asName)
	throws SQLException
	{
		sqlWhereSecHelper(prep, asName);
	}


	/******* GETTERS AND SETTERS ************/
	public String getMethylTablePrefix() {
		return methylTablePrefix;
	}


	public void setMethylTablePrefix(String methylTablePrefix) {
		this.methylTablePrefix = methylTablePrefix;
	}


	public String getFeatureTablePrefix() {
		return featureTablePrefix;
	}


	public void setFeatureTablePrefix(String featureTablePrefix) {
		this.featureTablePrefix = featureTablePrefix;
	}


	public int getMinCTreads() {
		return minCTreads;
	}


	public void setMinCTreads(int minCTreads) {
		this.minCTreads = minCTreads;
	}


	public boolean getUseNonconversionFilter() {
		return useNonconversionFilter;
	}


	public void setUseNonconversionFilter(boolean useNonconversionFilter) {
		this.useNonconversionFilter = useNonconversionFilter;
	}


	public double getMaxOppstrandAfrac() {
		return maxOppstrandAfrac;
	}


	public void setMaxOppstrandAfrac(double maxOppstrandAfrac) {
		this.maxOppstrandAfrac = maxOppstrandAfrac;
	}
	
	/******* BEHIND THE SCENES PRIVATE DB STUFF ********/

	/**
	 * @param prep If supplied, we not only return the SQL but fill in the prep with its correct values
	 * @return
	 */
	protected String sqlWhereSecHelper(PreparedStatement prep, String asName)
	throws SQLException
	{
		return this.sqlWhereSecHelper(prep, asName, 1);
	}
	
	protected String sqlWhereSecHelper(PreparedStatement prep, String asName, int curInd)
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
				String clause = String.format("(%schromPos>=? AND %schromPos<=?)", asSec, asSec);
				grClauses.add(clause);
				if (prep!=null)
				{
					prep.setInt(curInd++, gr.getStart());
					prep.setInt(curInd++, gr.getEnd());
				}
			}
			ListUtils.setDelim(" OR ");
			clauses.add("(" + ListUtils.excelLine(grClauses.toArray(new String[1])) + ")");
		}

		// feat filters
		if (this.featFilters.size()>0)
		{
//			for (MethylDbQuerier.FeatClass featFilter : this.featFilters)
//			{
//				// This is the best way to pick anything that OVERLAPS the range.
//				// Notice than CONTAINED WITHIN the range would be a different query
//				
//				// Add the flank to the SQL itself or as a param?  Seems like one script
//				// would use 1 or very few flank sizes, so i'll put it in SQL
//				
//				String clause = String.format("(!(%schromPosStart<? AND %schromPosEnd<?) AND !(%schromPosStart>? AND %schromPosEnd>?))", 
//						asSec, asSec,asSec,asSec);
//				clauses.add(clause);
//				if (prep!=null)
//				{
//					prep.setInt(curInd++, gr.getStart());
//					prep.setInt(curInd++, gr.getStart());
//					prep.setInt(curInd++, gr.getEnd());
//					prep.setInt(curInd++, gr.getEnd());
//				}
//			}
		}

		// minCTreads
		if (this.minCTreads > 0)
		{
			String cSec = (this.useNonconversionFilter) ? String.format("%scReads", asSec) :
				String.format("(%scReads+%scReadsNonconversionFilt)",asSec,asSec);
			String clause = String.format("((%s+%stReads) >= ?)",cSec, asSec);
			clauses.add(clause);
			if (prep != null)
			{
				prep.setInt(curInd++, this.minCTreads);
			}
		}
		
		
		// oppstrand gtoa ratio
		if (this.maxOppstrandAfrac < 1.0)
		{
			String clause = String.format("((%stotalReadsOpposite=0) OR ((%saReadsOpposite/%stotalReadsOpposite)<=?))",asSec, asSec, asSec);
			clauses.add(clause);
			if (prep != null)
			{
				prep.setDouble(curInd++, (double)this.maxOppstrandAfrac);
			}
		}
		
		ListUtils.setDelim(" AND ");
		String sql = ListUtils.excelLine(clauses.toArray(new String[1]));
		
		return sql;
	}
	
	public class FeatClass
	{
		String featType = null;
		int flank = 0;
		
		/**
		 * @param featType
		 * @param flank
		 */
		public FeatClass(String featType, int flank) {
			super();
			this.featType = featType;
			this.flank = flank;
		}
		
		
	}
	
	
}
