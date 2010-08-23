package edu.usc.epigenome.genomeLibs.MethylDb;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;

import org.biojava.bio.seq.StrandedFeature;
import org.biojava.bio.seq.StrandedFeature.Strand;

public class MethylRead implements Comparable<MethylRead> {

	public int readId = -1;
	public StrandedFeature.Strand strand = StrandedFeature.UNKNOWN;
	public TreeSet<Integer> methPositions = new TreeSet<Integer>();
	public TreeSet<Integer> unmethPositions = new TreeSet<Integer>();
	
	static StringBuffer c_sb = new StringBuffer(10000);

	
	/**
	 * @param strand
	 */
	public MethylRead(Strand inStrand, int inReadId) {
		super();
		this.strand = inStrand;
		this.readId = inReadId;
	}


	/****** OUTPUT *******/
	
	public static PrintWriter outputChromToFile(Set<MethylRead> cpgMap, String prefix, String sampleName, String chr, int minCphCoverage, double minCphMethFrac)
	throws IOException
	{
		
		String fn = prefix + sampleName + "_" + chr + ".txt";
		PrintWriter writer = new PrintWriter(new File(fn));
		
		outputMethylReadsToFile(writer, cpgMap, prefix, sampleName, chr, minCphCoverage, minCphMethFrac);
		
		//writer.close();
		return writer;
	}
	

	public static void outputMethylReadsToFile(PrintWriter pw, Set<MethylRead> cpgMap, String prefix, String sampleName, String chr, int minCphCoverage, double minCphMethFrac)
	throws IOException
	{
		//System.err.println("About to write " + cpgMap.size() + " Cytosines to file");
		Iterator<MethylRead> mrIt = cpgMap.iterator();
		MR: while (mrIt.hasNext())
		{
			MethylRead mr = mrIt.next();
			String line = mr.toString();
			pw.println(line);
		}
	}

	/****** info *******/
	
	public int numMeth()
	{
		return methPositions.size();
	}

	public int numUnmeth()
	{
		return unmethPositions.size();
	}
	
	public int numTotal()
	{
		return numMeth() + numUnmeth();
	}
	
	public double fracMeth()
	{
		return (double)numMeth()/(double)numTotal();
	}
	
	public int startPos()
	{
		int methMin = (this.numMeth()<1) ? Integer.MAX_VALUE : methPositions.first().intValue();
		int unmethMin = (this.numUnmeth()<1) ? Integer.MAX_VALUE : unmethPositions.first().intValue();
		
		return Math.min(methMin,unmethMin);
	}

	public int endPos()
	{
		int methMax = (this.numMeth()<1) ? Integer.MIN_VALUE : methPositions.last().intValue();
		int unmethMax = (this.numUnmeth()<1) ? Integer.MIN_VALUE : unmethPositions.last().intValue();
		
		return Math.max(methMax,unmethMax);
	}
	
	public int width()
	{
		return endPos() - startPos() + 1;
	}

	/****** EDITING *******/

	public void addPosition(int pos, boolean meth)
	{
		Set<Integer> set = (meth) ? methPositions : unmethPositions;
		set.add(new Integer(pos));
	}


	/****** COMPARISONS *******/

	
	public int compareTo(MethylRead o) {
		int out = 0;
		
		out = this.startPos()-o.startPos();
		if (out!=0) return out;
		
//		out = this.width()-o.width();
//		if (out!=0) return -out;

		out = (this.fracMeth()==o.fracMeth()) ? 0 : ( (this.fracMeth()>o.fracMeth()) ? 1 : -1);
		if (out!=0) return out;
		
		out = this.endPos()-o.endPos();
		if (out!=0) return out;
		
		out = this.readId-o.readId;
		
		return out;
	}



	@Override
	public boolean equals(Object obj) {
		if (MethylRead.class.isInstance(obj))
		{
			System.err.println("Using MethylRead.compareTo for MethylRead.equals()");
			return (this.compareTo((MethylRead)obj) == 0);
		}
		else
		{
			System.err.println("Using Object.equals for MethylRead.equals()");
			return super.equals(obj);
		}
	}


	/****** OUTPUT *******/
	
	/* (non-Javadoc)
	 * @see java.lang.Object#toString()
	 */
	@Override
	public String toString() {
		
		return String.format("Read (%c)\tnum=%d,\t%d-%d,\tmeth=%.0f\t[id=%d]",this.strand.getToken(),this.numTotal(),this.startPos(),this.endPos(),this.fracMeth()*100.0,this.readId);
		
	}
	
	


	
}
