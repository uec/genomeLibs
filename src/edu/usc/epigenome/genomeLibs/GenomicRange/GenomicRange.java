package edu.usc.epigenome.genomeLibs.GenomicRange;

import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.Iterator;
import java.util.Set;
import java.util.TreeSet;
import java.util.logging.Level;
import java.util.logging.Logger;

import org.biojava.bio.seq.StrandedFeature;
import org.usckeck.genome.ChromFeatures;

import edu.usc.epigenome.genomeLibs.BiojavaUtils;
import edu.usc.epigenome.genomeLibs.GoldAssembly;
import edu.usc.epigenome.genomeLibs.ChromStringComparator;
import edu.usc.epigenome.genomeLibs.ListUtils;


public class GenomicRange implements Cloneable, Comparable<GenomicRange> {
	
    
	//TODO Use org.biojavax.bio.seq.SimplePosition for position info
	
	/* Class vars */
	static final String NO_CHROM = "NO_CHROM";
	static final double NO_SCORE = Double.NEGATIVE_INFINITY;
	static final Comparator<GenomicRange> comparator = new GenomicRangeComparator(true);
	
	
	/* Obj vars */
	protected String chrom = NO_CHROM; // Should always be lower case (case-insensitive)
	protected int start = 0;
	protected int end = 0;
	protected StrandedFeature.Strand strand = StrandedFeature.UNKNOWN;
	protected double score = NO_SCORE;
	
	/* Constructors */
	public GenomicRange(String inChrom, int inStart, int inEnd)
	{
		this.setChrom(inChrom);
		this.setStart(inStart);
		this.setEnd(inEnd);
	}

	public GenomicRange(String inChrom, int inStart, int inEnd, StrandedFeature.Strand inStrand)
	{
		this.setChrom(inChrom);
		this.setStart(inStart);
		this.setEnd(inEnd);
		this.setStrand(inStrand);
	}


	// Make a non-coordinate, "fake" chrom that can be easily recognized.
	public static GenomicRange infiniteChromRange(String inChrom)
	{
		return new GenomicRange(inChrom, Integer.MIN_VALUE, Integer.MAX_VALUE);
	}
	
	public boolean isInfiniteChromRange()
	{
		return ((this.start==Integer.MIN_VALUE) && (this.end==Integer.MAX_VALUE));
	}


	
	public static GenomicRange fullChromRange(String inChrom, String genomeAssembly)
	{
		int end = Integer.MAX_VALUE;
		try
		{
			end = GoldAssembly.chromLengthStatic(inChrom, genomeAssembly);
		}
		catch (Exception e)
		{
			Logger.getLogger(Logger.GLOBAL_LOGGER_NAME).log(Level.WARNING,
					"Can't find chrom length for assembly " + genomeAssembly + ": " + e.toString());
		}
		return new GenomicRange(inChrom, 0, end);
	}
	
	/**
	 * @return the strand
	 */
	public StrandedFeature.Strand getStrand() {
		return strand;
	}




	/**
	 * @return the chrom
	 */
	public String getChrom() {
		
		String chr = this.chrom;

		chr = chr.replace('x', 'X');
		chr = chr.replace('y', 'Y');
		chr = chr.replace('m', 'M');
		
		return chr;
	}



	/**
	 * @param chrom the chrom to set
	 */
	public void setChrom(String chrom) {
		this.chrom = chrom.toLowerCase();
	}



	/**
	 * @return the start
	 */
	public int getStart() {
		return (int)start;
	}



	/**
	 * @param strand the strand to set
	 */
	public void setStrand(StrandedFeature.Strand strand) {
		this.strand = strand;
	}


	/**
	 * @param start the start to set
	 */
	public void setStart(int start) {
		this.start = start;
	}



	/**
	 * @return the end
	 */
	public int getEnd() {
		return (int)end;
	}



	/**
	 * @param end the end to set
	 */
	public void setEnd(int end) {
		this.end = end;
	}


	/**
	 * @return the score
	 */
	public double getScore() {
		return score;
	}

	public boolean overlaps(GenomicRange other)
	{
		// They have to be on the same chrom
		if (!this.getChrom().equalsIgnoreCase(other.getChrom())) return false;
		
		boolean overlap = !( (this.getStart()>other.getEnd()) || (other.getStart()>this.getEnd()) );
		return overlap;
	}

	/**
	 * @param score the score to set
	 */
	public void setScore(double score) {
		this.score = score;
	}

	/** Static util functions **/
	public static GenomicRange generateWindowFromCoord(String inChrom, int inCoord, int inWindSize)
	{
			return generateWindowFromCoord(inChrom, inCoord, inWindSize, StrandedFeature.UNKNOWN);
	}	
	
	public static GenomicRange generateWindowFromCoord(String inChrom, int inCoord, int inWindSize, StrandedFeature.Strand strand)
	{
		int s = ((int)Math.floor(inCoord/inWindSize))*inWindSize;
		int e = s + inWindSize - 1;
		return new GenomicRange(inChrom, s, e, strand);
	}
	
	/**
	 * @param inGenome
	 * @param inWindSize
	 * @param strandSpecific
	 * @param useSlowAndCorrectSorting
	 * @return
	 * @throws Exception
	 */
	public static Iterator<GenomicRange> allPossibleGenomicRanges(String inGenome, int inWindSize,
			boolean strandSpecific, boolean useSlowAndCorrectSorting)
	throws Exception
	{
		Set<GenomicRange> outSet = (useSlowAndCorrectSorting) ? 
				new TreeSet<GenomicRange>(new GenomicRangeComparator(false)) : new TreeSet<GenomicRange>();
		
		// Go through chroms
		Iterator<String> chromIt = GoldAssembly.chromIterator(inGenome, false);
		while (chromIt.hasNext())
		{
			String chrom = chromIt.next();
			int chrLen = GoldAssembly.chromLengthStatic(chrom, inGenome);
			
			for (int s = 0; s < chrLen; s+=inWindSize)
			{
				int e = s + inWindSize - 1;
				if (e > chrLen) e = chrLen;
				
				GenomicRange gr;
				if (strandSpecific)
				{
					gr = new GenomicRange(chrom, s, e, StrandedFeature.POSITIVE);
					outSet.add(gr);
					gr = new GenomicRange(chrom, s, e, StrandedFeature.NEGATIVE);
					outSet.add(gr);
				}					
				else
				{
					gr = new GenomicRange(chrom, s, e);
					outSet.add(gr);
				}
			}
		}

		System.err.println("Found " + outSet.size() + " genomic ranges size " + inWindSize);

		return outSet.iterator();
	}
	

	/** Cloneable **/
	
	public GenomicRange clone()
	{
		GenomicRange newInst = null;
		try
		{
			newInst = (GenomicRange)super.clone();
		}
		catch (CloneNotSupportedException e)
		{
			// Should never get here
			e.printStackTrace();
			System.exit(0);
		}
		
		return newInst;
	}

	/** 
	 * If symbol and strand are equivalent, the two are equal
	 */
	@Override
	public boolean equals(Object obj) {
		if (this == obj) return true;
		
		if (obj instanceof GenomicRange)
		{
			return (this.compareTo((GenomicRange)obj) == 0);
		}
		else
		{	
			return false;
		}
	}

	/**
	 * If symbol and strand are equivalent, the two are equal
	 */
	@Override
	public int hashCode() {
		return this.equalsKey().hashCode();
	}
	
	protected String equalsKey()
	{
		return commaSeparatedLine();
	}
		
	/* (non-Javadoc)
	 * @see java.lang.Comparable#compareTo(java.lang.Object)
	 */
	public int compareTo(GenomicRange o) {
		return comparator.compare(this, o);
	}

	public String commaSeparatedLine()
	{
		return this.commaSeparatedLine(false);
	}	

	public String gffStr()
	{

		String[] cols = {
					this.getChrom(),
					"DEFAULT",
					"exon",
					String.format("%d",this.getStart()),
					String.format("%d",this.getEnd()),
					String.format("%.1f",this.getScore()),
					Character.toString(BiojavaUtils.strandToSymbol(this.getStrand())),
					"."
		};
		
		return ListUtils.tabbedLine(Arrays.asList(cols));
		
	}	

	public String commaSeparatedLine(boolean intChrom)
	{
		String chr = this.getChrom();
		if (intChrom)
		{
			chr = Integer.toString((new ChromFeatures()).chrom_from_public_str(chr));
		}
		return String.format("%s,%d,%d,%d",chr, BiojavaUtils.strandToInt(this.getStrand()), 
				this.getStart(), this.getEnd());
	}

	/* (non-Javadoc)
	 * @see java.lang.Object#toString()
	 */
	@Override
	public String toString() {
		return this.commaSeparatedLine();
	}
	

	
	
}
