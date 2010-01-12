package edu.usc.epigenome.genomeLibs.MethylDb.CpgSummarizers;

import java.util.List;
import java.util.Set;
import java.util.logging.Logger;

import com.googlecode.charts4j.Color;

import edu.usc.epigenome.genomeLibs.FeatDb.FeatDbQuerier;
import edu.usc.epigenome.genomeLibs.FeatDb.FeatIterator;
import edu.usc.epigenome.genomeLibs.GenomicRange.GenomicRange;
import edu.usc.epigenome.genomeLibs.MethylDb.Cpg;
import edu.usc.epigenome.genomeLibs.MethylDb.MethylDbQuerier;

public class CpgDensitySummarizer extends CpgSummarizer {

	
	
	
	/**
	 * 
	 */
	public CpgDensitySummarizer() {
		super();
		// TODO Auto-generated constructor stub
	}

	/**
	 * @param inQuerier
	 */
	public CpgDensitySummarizer(MethylDbQuerier inQuerier) {
		super(inQuerier);
		// TODO Auto-generated constructor stub
	}

	@Override
	public void init(MethylDbQuerier inQuerier) {
		// TODO Auto-generated method stub
		super.init(inQuerier);
		
		this.setColor(Color.RED);
		this.setNaturalMax(1.0);
		this.setNaturalMin(0.0);
		
		this.setLabel("%Cpg");
	

		this.valsMin = 0.0;
		this.valsMax = 1.0;
		this.numVals = 0.0;
	}
	
	// Must call this every time you set up new genomic ranges.
	public void addRanges(MethylDbQuerier inQuerier)
	{
		/*********** TO DO
		 * TO DO
		 * Check if we want union or intersection
		 */
		try
		{
			List<String> featTypes = inQuerier.getFeatFilterTypes();
			Set<GenomicRange> ranges = inQuerier.getRangeFilters();
			this.numVals += FeatIterator.TotalFeatSize(featTypes, ranges, true);
			Logger.getLogger(Logger.GLOBAL_LOGGER_NAME).severe("Total feat size for " + inQuerier.getFeatFilterTypesList() + " = " + this.numVals);
		}
		catch (Exception e)
		{
			Logger.getLogger(Logger.GLOBAL_LOGGER_NAME).severe("Can't get total feat size:\n"+e.getMessage());
			e.printStackTrace();
		}
	}

		
	@Override
	public void streamCpg(Cpg cpg) {
		// TODO Auto-generated method stub
		super.streamCpg(cpg);

		// Count every CpG
		this.valsTotal++;
		this.valsSquareTotal++;
	}


}
