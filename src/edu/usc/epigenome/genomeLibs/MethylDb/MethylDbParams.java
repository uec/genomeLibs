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

public class MethylDbParams {

	public final static String DEFAULT_TABLE_PREFIX = "methylCGsRich_tumor10x_";
	public final static String DEFAULT_CONN_STR = "jdbc:mysql://localhost/cr?user=benb";
	
	
	public String tablePrefix = DEFAULT_TABLE_PREFIX;
	public String connStr = DEFAULT_CONN_STR;
	
	
	// Ranges
	public Set<GenomicRange> rangeFilters = new HashSet<GenomicRange>(1);
	public Set<String> featFilters = new HashSet<String>(1);

	
	// Filtering
	public int minCTreads = 0;
	public boolean useNonconversionFilter = true;
	public double maxOppstrandAfrac = 0.2; // Double.MAX_VALUE;
	
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


	
	/******* DB STUFF ********/
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

	
	/**
	 * @param prep If supplied, we not only return the SQL but fill in the prep with its correct values
	 * @return
	 */
	public String sqlWhereSecHelper(PreparedStatement prep, String asName)
	throws SQLException
	{
		return this.sqlWhereSecHelper(prep, asName, 1);
	}
	
	public String sqlWhereSecHelper(PreparedStatement prep, String asName, int curInd)
	throws SQLException
	{
		String asSec = (asName==null) ? "" : (asName+".");

		List<String> clauses = new ArrayList<String>(20);
		
		
		// Range filters
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
		
		// feat filters
		for (String featTablePrefix : this.featFilters)
		{
			
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
	
}
