package edu.usc.epigenome.genomeLibs.MethylDb;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Iterator;
import java.util.Map;

import org.biojava.bio.seq.StrandedFeature;

public class Cpg implements Comparable {

	// Object vars
	public int chromPos = 0;
	public boolean negStrand = false;
	public short totalReads = 0;
	public short cReads = 0;
	public short cReadsNonconversionFilt = 0;
	public short tReads = 0;
	public short agReads = 0;
	public short totalReadsOpposite = 0;
	public short aReadsOpposite = 0;



	/**
	 * @param chromPos
	 * @param negStrand
	 * @param totalReads
	 * @param cReads
	 * @param cReadsNonconversionFilt
	 * @param tReads
	 * @param agReads
	 * @param totalReadsOpposite
	 * @param aReadsOpposite
	 */
	public Cpg(int chromPos, boolean negStrand) {
		super();
		this.chromPos = chromPos;
		this.negStrand = negStrand;
	}
	
	public Cpg(int chromPos, boolean negStrand, short totalReads, short cReads,
			short cReadsNonconversionFilt, short tReads, short agReads,
			short totalReadsOpposite, short aReadsOpposite) {
		super();
		this.chromPos = chromPos;
		this.negStrand = negStrand;
		this.totalReads = totalReads;
		this.cReads = cReads;
		this.cReadsNonconversionFilt = cReadsNonconversionFilt;
		this.tReads = tReads;
		this.agReads = agReads;
		this.totalReadsOpposite = totalReadsOpposite;
		this.aReadsOpposite = aReadsOpposite;
	}
	
	/*** Overridden Comparable methods
	 */

	/* (non-Javadoc)
	 * @see java.lang.Comparable#compareTo(java.lang.Object)
	 */
	public int compareTo(Object o) {

		Cpg other = (Cpg)o;
		
		return (this.chromPos - other.chromPos);
	}


	/* (non-Javadoc)
	 * @see java.lang.Object#equals(java.lang.Object)
	 */
	@Override
	public boolean equals(Object obj) {
		return 	(this.compareTo(obj) == 0);
	}
	
	
	public double fracMeth(boolean useNonconvFilt)
	{
		int cs = (useNonconvFilt) ? this.cReads : (this.cReads + this.cReadsNonconversionFilt); 
		return (double)cs / ((double)cs + (double)this.tReads);
	}

	public String variableStepWigLine(boolean useNonconvFilter)
	{
		return String.format("%d\t%f",chromPos, this.fracMeth(useNonconvFilter));
	}

	public static void outputChromToFile(Map<Integer,Cpg> cpgMap, String prefix, String sampleName, String chr)
	throws IOException
	{
		
		String fn = prefix + sampleName + "_" + chr + ".txt";
		PrintWriter writer = new PrintWriter(new File(fn));
		
		Iterator<Cpg> cpgIt = cpgMap.values().iterator();
		while (cpgIt.hasNext())
		{
			Cpg cpg = cpgIt.next();
			String line = cpg.toString();
			writer.println(line);
		}
		writer.close();
	}

	@Override
	public String toString() {

		return String.format("%d\t%c\t%d\t%d\t%d\t%d\t%d\t%d\t%d", 
				chromPos,
				(negStrand) ? '-' : '+',
				totalReads,
				cReads,
				cReadsNonconversionFilt,
				tReads,
				agReads,
				totalReadsOpposite,
				aReadsOpposite);
	}
	
	public StrandedFeature.Strand getStrand()
	{
		return (this.negStrand) ? StrandedFeature.NEGATIVE : StrandedFeature.POSITIVE;
	}
	
	public String getStrandStr()
	{
		return (this.negStrand) ? "-" : "+";
	}
	
	public short totalReadsC(boolean useNonconvFilt)
	{
		int cs = (useNonconvFilt) ? this.cReads : (this.cReads + this.cReadsNonconversionFilt); 
		return (short)cs;
	}
	
	public short totalReadsCorT(boolean useNonconvFilt)
	{
		short cs = totalReadsC(useNonconvFilt);
		return (short)(cs + this.tReads);
	}
	
	
	
	
}
