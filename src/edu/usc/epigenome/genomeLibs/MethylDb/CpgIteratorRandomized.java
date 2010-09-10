package edu.usc.epigenome.genomeLibs.MethylDb;

import java.util.ArrayList;
import java.util.Collections;

/**
 * @author benb
 *
 * Outputs all the CpGs of the underlying CpG iterator, but scrambles the methylation values 
 * 
 */
public class CpgIteratorRandomized extends CpgIterator {

	
	protected ArrayList<Cpg> cpgsRandomized = null;
	protected int listPos = -1;
	
	public CpgIteratorRandomized(MethylDbQuerier inParams)
	throws Exception 
	{
		this.init(inParams);
	}



	public CpgIteratorRandomized(String chrom, int startCoord, int endCoord, String inTablePrefix)
	throws Exception 
	{
		params = new MethylDbQuerier();
		if (inTablePrefix!=null) params.methylTablePrefix = inTablePrefix;
		params.addRangeFilter(chrom, startCoord, endCoord);
		this.init(params);
	}
		
	
	public CpgIteratorRandomized(String chrom, String inTablePrefix)
	throws Exception 
	{
		params = new MethylDbQuerier();
		if (inTablePrefix!=null) params.methylTablePrefix = inTablePrefix;
		params.addRangeFilter(chrom);
		this.init(params);

	}

	public CpgIteratorRandomized(MethylDbQuerier inParams, String inTablePrefix)
	throws Exception 
	{
		params = inParams;
		if (inTablePrefix!=null) params.methylTablePrefix = inTablePrefix;
		this.init(params);
	}	
	


	/* (non-Javadoc)
	 * @see edu.usc.epigenome.genomeLibs.MethylDb.CpgIterator#init(edu.usc.epigenome.genomeLibs.MethylDb.MethylDbQuerier)
	 */
	@Override
	public int init(MethylDbQuerier inParams) throws Exception {
		
		// This is where we run through and fill in our
		int numRows = super.init(inParams);
		
		System.err.printf("(%s) About to load %d meth values\n", this, numRows);
		this.cpgsRandomized = new ArrayList<Cpg>(numRows);
		int i = 1;
		while (super.hasNext())
		{
			Cpg cpg = super.next();
			//this.methValsRandomized.add(cpg.fracMeth(params.useNonconversionFilter));
			this.cpgsRandomized.add(cpg);
			if ((i%1E5)==0) System.err.printf("\t\tOn CpG %d\n",i);
		}
		System.err.printf("\tLoaded %d meth values, about to permute\n", numRows);
		Collections.shuffle(this.cpgsRandomized);
		System.err.printf("\tFinished permuting list of %d meth values\n", numRows);
		
		// Initialize list position
		this.listPos = 0;
		
		// Now we run through a second time, to be passed on as the iterator. 
		return super.init(inParams);
	}

	/* (non-Javadoc)
	 * @see edu.usc.epigenome.genomeLibs.MethylDb.CpgIterator#next()
	 */
	@Override
	public Cpg next() {
		// Return the randomized version, except use our actual chrom position and strand
		Cpg actual = super.next();
		Cpg randomized = this.cpgsRandomized.get(this.listPos++);
//		if (actual.fracMeth(true)<0.05)
//			System.err.printf("About to get randomized CpG #%d: Actual=%.2f\tRandomized=%.2f\n", this.listPos, actual.fracMeth(true),randomized.fracMeth(true));
		
		randomized.chromPos = actual.chromPos;
		randomized.negStrand = actual.negStrand;
		randomized.cpgWeight = actual.cpgWeight;
		
		return randomized;
	}

	/* (non-Javadoc)
	 * @see edu.usc.epigenome.genomeLibs.MethylDb.CpgIterator#remove()
	 */
	@Override
	public void remove() {
		// TODO Auto-generated method stub
		super.remove();
	}
	
	
	

}
