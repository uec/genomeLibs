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
import edu.usc.epigenome.genomeLibs.FeatDb.FeatDbQuerier;
import edu.usc.epigenome.genomeLibs.GenomicRange.GenomicRange;

public class MethylDbQuerier {


	public final static String DEFAULT_METHYL_TABLE_PREFIX = "methylCGsRich_tumor10x_";
	public final static String DEFAULT_FEAT_TABLE_PREFIX = "features_";
	public final static String DEFAULT_CONN_STR = "jdbc:mysql://localhost/cr?user=benb";

	public final static String connStr = DEFAULT_CONN_STR;
	
	
	public String methylTablePrefix = DEFAULT_METHYL_TABLE_PREFIX;
	public String featureTablePrefix = DEFAULT_FEAT_TABLE_PREFIX;
	
	
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
	
	public void clearRangeFilters()
	{
		this.rangeFilters = new HashSet<GenomicRange>(1);
	}
	
	public void addRangeFilter(String chrom)
	{
		// Set this up in a way we can detect this later, because we will want
		// to remove it for database efficiency.
		GenomicRange gr = GenomicRange.infiniteChromRange(chrom);
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
	
	public String getMethylTable()
	throws Exception
	{
		String chrom = this.getChrom();
		String table = this.methylTablePrefix + chrom;

		return table;
	}

	public String getFeatTable()
	throws Exception
	{
		String chrom = this.getChrom();
		String table = this.featureTablePrefix + chrom;

		return table;
	}

	
	
	
	



	/******* GETTERS AND SETTERS ************/
	
	public boolean usesFeatTable()
	{
		return (this.featFilters.size() > 0);
	}
	
	public List<String> getFeatTableNames()
	{
		//int nF = this.featFilters.size();
		List<String> out = new ArrayList<String>();
		for (MethylDbQuerier.FeatClass feat : this.featFilters)
		{
			out.add(feat.featType);
		}

		return out;
	}
	
	public String getFeatTableNameList()
	{
		List<String> tables = this.getFeatTableNames();
		ListUtils.setDelim(", ");
		return ListUtils.excelLine(tables);
	}
	
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
			MethylDbQuerier.HelperOutput output = sqlWhereSecHelper(null, asName);
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
		MethylDbQuerier.HelperOutput output = sqlWhereSecHelper(prep, asName);
		return output.newCurInd;
	}
	
	
	/******* BEHIND THE SCENES PRIVATE DB STUFF ********/

	/**
	 * @param prep If supplied, we not only return the SQL but fill in the prep with its correct values
	 * @return
	 */
	protected MethylDbQuerier.HelperOutput sqlWhereSecHelper(PreparedStatement prep, String asName)
	throws SQLException
	{
		return this.sqlWhereSecHelper(prep, asName, 1);
	}
	
	protected MethylDbQuerier.HelperOutput sqlWhereSecHelper(PreparedStatement prep, String asName, int curInd)
	throws SQLException
	{
		String asSec = (asName==null) ? "" : (asName+".");

		List<String> clauses = new ArrayList<String>(20);
		
		
		// Range filters
		if (this.rangeFilters.size()>0)
		{
			// Testing with the indices in my table shows that if you're doing a whole chromosome, you 
			// NEVER want to put artificial filters on (like the "whole chromosome" filter we sometimes
			// use.  So let's detect this.
//			mysql> select straight_join count(*) from features_chr11, methylCGsRich_normal123009_chr11 cpg  WHERE ((cpg.chromPos>=0 AND cpg.chromPos<=5531960706)) AND  ((cpg.chromPos>=(chromPosStart)) AND (cpg.chromPos<=(chromPosEnd))  AND (featType = 'hg18.ES.H3K27me3.HMM.gtf'))  AND ((cpg.cReads+cpg.tReads) >= 4) AND ((cpg.totalReadsOpposite=0) OR ((cpg.aReadsOpposite/cpg.totalReadsOpposite)<=0.2)) AND (chromPosEnd>0) AND (chromPosStart<5531960706) ORDER BY chromPos ;
//			+----------+
//			| count(*) |
//			+----------+
//			|    49327 |
//			+----------+
//			1 row in set (1 min 51.58 sec)
//
//			mysql> select straight_join count(*) from features_chr11, methylCGsRich_normal123009_chr11 cpg  WHERE  ((cpg.chromPos>=(chromPosStart)) AND (cpg.chromPos<=(chromPosEnd))  AND (featType = 'hg18.ES.H3K27me3.HMM.gtf'))  AND ((cpg.cReads+cpg.tReads) >= 4) AND ((cpg.totalReadsOpposite=0) OR ((cpg.aReadsOpposite/cpg.totalReadsOpposite)<=0.2)) ORDER BY chromPos ;
//			+----------+
//			| count(*) |
//			+----------+
//			|    49327 |
//			+----------+
//			1 row in set (0.11 sec)
			
			List<String> grClauses = new ArrayList<String>(this.rangeFilters.size());
			for (GenomicRange gr : this.rangeFilters)
			{
				// Rule out infinite ones
				if (!gr.isInfiniteChromRange())
				{
					String clause = String.format("(%schromPos>=? AND %schromPos<=?)", asSec, asSec);
					grClauses.add(clause);
					if (prep!=null)
					{
						prep.setInt(curInd++, gr.getStart());
						prep.setInt(curInd++, gr.getEnd());
					}
				}
			}
			
			if (grClauses.size()>0)
			{
				ListUtils.setDelim(" OR ");
				clauses.add("(" + ListUtils.excelLine(grClauses.toArray(new String[1])) + ")");
			}
		}

		// feat filters
		if (this.featFilters.size()>0)
		{
			for (MethylDbQuerier.FeatClass featFilter : this.featFilters)
			{
				// When we have a feature filter, we generally use that index first since the 
				// number of rows is much smaller in the feature table than the meth table.
				// So if we add a redundant coordinate filter for the features, we can often
				// speed up the query incredibly:

//				mysql> select straight_join count(*) from features_chr11, methylCGsRich_normal123009_chr11 cpg  WHERE ((cpg.chromPos>=30683099 AND cpg.chromPos<=31960706)) AND  ((cpg.chromPos>=(chromPosStart)) AND (cpg.chromPos<=(chromPosEnd))  AND (featType = 'hg18.ES.H3K27me3.HMM.gtf'))  AND ((cpg.cReads+cpg.tReads) >= 4) AND ((cpg.totalReadsOpposite=0) OR ((cpg.aReadsOpposite/cpg.totalReadsOpposite)<=0.2))  ORDER BY chromPos ;
//				+----------+
//				| count(*) |
//				+----------+
//				|     2210 |
//				+----------+
//				1 row in set (0.82 sec)
//
//				mysql> select straight_join count(*) from features_chr11, methylCGsRich_normal123009_chr11 cpg  WHERE ((cpg.chromPos>=30683099 AND cpg.chromPos<=31960706)) AND  ((cpg.chromPos>=(chromPosStart)) AND (cpg.chromPos<=(chromPosEnd))  AND (featType = 'hg18.ES.H3K27me3.HMM.gtf'))  AND ((cpg.cReads+cpg.tReads) >= 4) AND ((cpg.totalReadsOpposite=0) OR ((cpg.aReadsOpposite/cpg.totalReadsOpposite)<=0.2)) AND (chromPosEnd>30683099) AND (chromPosStart<31960706) ORDER BY chromPos ;
//				+----------+
//				| count(*) |
//				+----------+
//				|     2210 |
//				+----------+
//				1 row in set (0.02 sec)

				
				// Since this is the case, the easiest thing to do is to add in the section from the feat DB helper.
				FeatDbQuerier featQuery = new FeatDbQuerier();
				featQuery.addFeatFilter(featFilter.featType);
				for (GenomicRange gr : this.rangeFilters)
				{
					if (!gr.isInfiniteChromRange()) featQuery.addRangeFilter(gr);
				}
				FeatDbQuerier.HelperOutput featOutput = featQuery.sqlWhereSecHelper(prep, null, curInd);
				//System.err.printf("Adding featOutput to query.  curInd before:%d, curInd after:%d\n",curInd,featOutput.newCurInd);
				curInd = featOutput.newCurInd;
				clauses.add(featOutput.sql);
				
				// Then join to meth table.
				String clause = String.format(" ((%schromPos>=(chromPosStart-%d)) AND (%schromPos<=(chromPosEnd+%d))) ",
						asSec, featFilter.flank, asSec, featFilter.flank);
				clauses.add(clause);
			}
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
		
		return new MethylDbQuerier.HelperOutput(sql, curInd);
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
