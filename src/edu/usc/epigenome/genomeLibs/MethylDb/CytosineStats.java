package edu.usc.epigenome.genomeLibs.MethylDb;

import java.util.ArrayList;
import java.util.List;

import edu.usc.epigenome.genomeLibs.MatUtils;

public class CytosineStats {

	private static final int MAXCOVERAGE = 100000;
	
	protected double maxOppstrandAfrac = 0.0;
	protected double maxNextNonGfrac = 0.0;

	// [0] = CpG, [1] = CpH, [2] = cytosine, uncertain context, [3] potential non-cytosine
	private static final int CAT_CG = 0;
	private static final int CAT_CH = 1;
	private static final int CAT_CYT_UNCERTAIN = 2;
	private static final int CAT_POTENTIAL_NONCYT = 3;
	
	
	protected int uniqueCountsByCvg[][] = null;
	protected int totalReads[] = null;
	protected int totalCreads[] = null;
	protected int totalTreads[] = null;
	protected double totalPercentC[] = null;
	protected int totalPercentCcount[] = null;
	
	/**
	 * @param maxOppstrandAfrac
	 * @param maxNextNonGfrac
	 */
	public CytosineStats(double maxOppstrandAfrac, double maxNextNonGfrac) {
		super();
		this.maxOppstrandAfrac = maxOppstrandAfrac;
		this.maxNextNonGfrac = maxNextNonGfrac;
		initArrays();
	}
	
	protected void initArrays()
	{
		uniqueCountsByCvg = new int[4][MAXCOVERAGE];
		totalReads = new int[4];
		totalCreads = new int[4];
		totalTreads = new int[4];
		totalPercentC = new double[4];
		totalPercentCcount = new int[4];
	}
	
	public void streamCytosine(Cpg cyt)
	{
		// Categorize
		int cat = -1;
		if (cyt.totalReads==0)
		{
			// Some have reads only on the opposite strand, but no coverage
			return;
		}
		else if (cyt.fracOppositeA() > this.maxOppstrandAfrac)
		{
			cat = CAT_POTENTIAL_NONCYT;
		}
		else if (cyt.fracReadsCorT() < (1-this.maxOppstrandAfrac))
		{
			cat = CAT_POTENTIAL_NONCYT;
		}
		else if (cyt.fracNextBaseG() >= (1-this.maxNextNonGfrac))
		{
			cat = CAT_CG;
		}
		else if (cyt.fracNextBaseG() < this.maxNextNonGfrac)
		{
			cat = CAT_CH;	
		}
		else
		{
			cat = CAT_CYT_UNCERTAIN;
		}
		
		if (cyt.totalReads == 0) 
		{
			System.err.printf("Following has 0 reads: %s\n",cyt.toString());
		}
		
		uniqueCountsByCvg[cat][cyt.totalReads]++;
		totalReads[cat] += cyt.totalReads-cyt.cReadsNonconversionFilt;
		totalCreads[cat] += cyt.cReadsNonconversionFilt;
		totalTreads[cat] += cyt.tReads;
		
		double fracC = cyt.fracMeth(true);
		if (!Double.isNaN(fracC))
		{
			totalPercentC[cat] += fracC;
			totalPercentCcount[cat]++;
		}
	
	}
	
	public List<String> headerStrings()
	{
		List<String> out = new ArrayList<String>();
		
		for (int cat = 0; cat < 4; cat++)
		{
			out.add(String.format("unique-%s", catToStr(cat)));
			out.add(String.format("totalReads-%s", catToStr(cat)));
			out.add(String.format("gte4-%s", catToStr(cat)));
			out.add(String.format("gte5-%s", catToStr(cat)));
			out.add(String.format("gte10-%s", catToStr(cat)));
			out.add(String.format("meanPercentT-%s", catToStr(cat)));
		}
	
		for (int cat = 0; cat < 4; cat++)
		{
			out.add(String.format("totalCTreads-%s", catToStr(cat)));
			out.add(String.format("totalCreads-%s", catToStr(cat)));
			out.add(String.format("e0-%s", catToStr(cat)));
			out.add(String.format("e1-%s", catToStr(cat)));
			out.add(String.format("e2-%s", catToStr(cat)));
		}
		
		return out;
	}
	
	public List<String> summaryStrings()
	{
		List<String> out = new ArrayList<String>();
		for (int cat = 0; cat < 4; cat++)
		{
			int totalUnique = MatUtils.nanSum(uniqueCountsByCvg[cat]);
			out.add(String.format("%d", totalUnique));
			out.add(String.format("%d", totalReads[cat]));
			out.add(String.format("%d", MatUtils.nanSum(uniqueCountsByCvg[cat],4,uniqueCountsByCvg[cat].length-4)));
			out.add(String.format("%d", MatUtils.nanSum(uniqueCountsByCvg[cat],5,uniqueCountsByCvg[cat].length-5)));
			out.add(String.format("%d", MatUtils.nanSum(uniqueCountsByCvg[cat],10,uniqueCountsByCvg[cat].length-10)));
			out.add(String.format("%.4f", 1.0-(totalPercentC[cat]/(double)totalPercentCcount[cat])));
		}
		for (int cat = 0; cat < 4; cat++)
		{
			out.add(String.format("%d", totalCreads[cat]+totalTreads[cat]));
			out.add(String.format("%d", totalCreads[cat]));
			out.add(String.format("%d", uniqueCountsByCvg[cat][0]));
			out.add(String.format("%d", uniqueCountsByCvg[cat][1]));
			out.add(String.format("%d", uniqueCountsByCvg[cat][2]));
		}
		
		return out;
	}
	
	protected String catToStr(int cat)
	{
		String out = null;
		switch (cat)
		{
		case CAT_POTENTIAL_NONCYT: out="possibleNonCytosine"; break;
		case CAT_CYT_UNCERTAIN: out="cytosineUnknownCG"; break;
		case CAT_CG: out="CG"; break;
		case CAT_CH: out="CH"; break;
		default: System.err.printf("catToStr only accepts vals from 0-3, got %d\n", cat); System.exit(1); break;
		}
		return out;
	}
	
}
