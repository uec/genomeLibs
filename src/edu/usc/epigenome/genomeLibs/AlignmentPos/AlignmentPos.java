package edu.usc.epigenome.genomeLibs.AlignmentPos;

import java.util.*;

import edu.usc.epigenome.genomeLibs.GFFUtils;
import edu.usc.epigenome.genomeLibs.Counters.SymbolCounter;
import edu.usc.epigenome.genomeLibs.Counters.SymbolCounterStratified;
import edu.usc.epigenome.genomeLibs.ReadPos.ReadPos;

import org.biojava.bio.program.gff.*;
import org.biojava.bio.seq.StrandedFeature;
import org.biojava.bio.symbol.*;
import org.biojava.bio.seq.DNATools;

public abstract class AlignmentPos implements Cloneable {

	/* Class vars */
	public static final double NO_COVERAGE = -1.0;
	public static final String NULL_CHROM = "NULLCHROM";
	public static final AlignmentPosOptions DEFAULT_AP_OPTIONS = new AlignmentPosOptions();

	/* Obj vars */
	protected Symbol ref;
	protected String chr;
	protected int pos;
	protected StrandedFeature.Strand strand = StrandedFeature.UNKNOWN;
	protected AlignmentPosOptions apOptions = null;
	


	/*****************
	 *  Constructors
	 */
	public AlignmentPos()
	{
	}
		
	public AlignmentPos(char inRef, String inChr, int inPos, AlignmentPosOptions inApOptions)
	{
		this.setRef(inRef);
		chr = inChr;
		pos = inPos;
		apOptions = inApOptions;
	}

	public AlignmentPos(Symbol inRef, String inChr, int inPos, AlignmentPosOptions inApOptions)
	{
		ref = inRef;
		chr = inChr;
		pos = inPos;
		apOptions = inApOptions;
	}


	public AlignmentPos(AlignmentPos ap)
	{
		ref = ap.ref;
		chr = ap.chr;
		pos = ap.pos;
		apOptions = ap.apOptions;
		strand = ap.strand;
	}


	/**
	 * 
	 * Accessors
	 * 
	 */

	/**
	 * @return the chr
	 */
	public String getChr() {
		return chr;
	}

	/**
	 * @param chr the chr to set
	 */
	public void setChr(String chr) {
		this.chr = chr;
	}

	/**
	 * @return the apOptions
	 */
	public AlignmentPosOptions getApOptions() {
		return apOptions;
	}

	/**
	 * @param apOptions the apOptions to set
	 */
	public void setApOptions(AlignmentPosOptions inApOptions) {
		this.apOptions = inApOptions;
	}

	/**
	 * @param ref the ref to set
	 */
	public void setRef(Symbol inRef) {
		this.ref = inRef;
	}

	/**
	 * @param ref the ref to set
	 */
	public void setRef(char inRef) {
		try
		{
			this.ref = DNATools.forSymbol(inRef);
		}
		catch (IllegalSymbolException e)
		{
			System.err.println(e);
		}
	}
	
	/**
	 * @return the pos
	 */
	public int getPos() {
		return pos;
	}

	/**
	 * @param pos the pos to set
	 */
	public void setPos(int pos) {
		this.pos = pos;
	}

	/**
	 * @return the strand
	 */
	public StrandedFeature.Strand getStrand() {
		return strand;
	}

	/**
	 * @return the strand symbol [+-.]
	 */
	public char getStrandSymbol() {
		StrandedFeature.Strand s = this.getStrand();
		
		if (s == StrandedFeature.POSITIVE) return '+';
		if (s == StrandedFeature.NEGATIVE) return '-';
		
		return '.';
	}

	/**
	 * @return the strand integer (-1, 0, +1)
	 */
	public int getStrandInt() {
		StrandedFeature.Strand s = this.getStrand();
		
		if (s == StrandedFeature.POSITIVE) return 1;
		if (s == StrandedFeature.NEGATIVE) return -1;
		
		return 0;
	}

	/**
	 * @param strand the strand to set
	 */
	public void setStrand(StrandedFeature.Strand strand) {
		this.strand = strand;
	}

	public static String getRefTokens(Collection<AlignmentPos> aps)
	{
		int len = aps.size();
		String out = "";
		if (len > 0)
		{
			AlignmentPos[] ar = new AlignmentPos[len];
			aps.toArray(ar);
			out = getRefTokens(ar); 
		}
		
		return out; 
	}
		
	public static String getRefTokens(AlignmentPos[] aps)
	{
		StringBuffer buf = new StringBuffer(aps.length);
		for (int i = 0; i < aps.length; i++)
		{
			buf.append(aps[i].getRefToken());
		}
		
		return buf.toString();
	}

	public static Symbol[] getSymbols(AlignmentPos[] aps)
	{
		
		Symbol[] out= new Symbol[aps.length];
		for (int i = 0; i < aps.length; i++)
		{
			out[i] = aps[i].getRef();
		}
		
		return out;
	}

	
	public static String getPosString(Collection<AlignmentPos> aps)
	{
		int len = aps.size();
		String out = "";
		if (len > 0)
		{
			AlignmentPos[] ar = new AlignmentPos[len];
			aps.toArray(ar);
			out = getPosString(ar); 
		}
		
		return out; 
	}

	static String getPosString(AlignmentPos[] aps)
	{
		StringBuffer buf = new StringBuffer(aps.length*3);
		for (int i = 0; i < aps.length; i++)
		{
			buf.append( aps[i].getPos() + ",");
		}
		
		return buf.toString();
	}

	public static String getCpgDensityStr(AlignmentPos[] aps)
	{
		double dens = getCpgDensity(aps);
		return String.format("%.2f", dens);
	}
	
	public static double getCpgDensity(AlignmentPos[] aps)
	{
		int totalCpgs = 0;
		int total = 0;
		for (int i = 1; i < aps.length; i++)
		{
			Symbol pre = aps[i-1].getRef();
			Symbol post = aps[i].getRef();
			
			if (pre != DNATools.n())
			{
				total++;
				if ((pre == DNATools.c()) && (post == DNATools.g())) totalCpgs++;
			}
		}
		return (double)totalCpgs/(double)total; 
	}
	
	public char getRefToken()
	{
		return this.getRefToken(true);
	}
	
	public char getRefToken(boolean reference_forward_strand)
	{
		Symbol ref = this.getRef(reference_forward_strand);
		char c = '0';
		try
		{
			c = DNATools.dnaToken(ref);
		}
		catch (Exception e)
		{
			e.printStackTrace();
			System.exit(0);
		}
		return c;
	}

	public Symbol getRef(boolean reference_forward_strand)
	{
		Symbol c = this.ref;
		try
		{
			if (!reference_forward_strand) c = DNATools.complement(c);
		}
		catch (Exception e)
		{
			e.printStackTrace();
			System.exit(0);
		}
		return c;
	}

	public Symbol getRef()
	{
		return this.getRef(true);
	}

	protected Symbol getRefFlipped()
	{
		return getRef(false);
	}
	
	public int getTotalDepth()
	{
		int[] d = getDepth();
		return d[0] + d[1];
	}
	
	public int getDepth(boolean fw)
	{
		return getDepth(true, fw);
	}
	
	public int getDepthNoIdentical(boolean fw)
	{
		int orig = this.apOptions.maxIdentical;
		this.apOptions.maxIdentical = 1;
		int depth = getDepth(true, fw);
		this.apOptions.maxIdentical = orig;
		return depth;
	}

	public int getDepthWithIdentical(boolean fw)
	{
		int orig = this.apOptions.maxIdentical;
		this.apOptions.maxIdentical = 0;
		int depth = getDepth(true, fw);
		this.apOptions.maxIdentical = orig;
		return depth;
	}
	
	public int getDepth(boolean reference_forward_strand, boolean read_same_strand)
	{
		int[] d = getDepth();
		return d[ (reference_forward_strand == read_same_strand) ? 0 : 1 ];
	}
	
	
	
	abstract public int[] getDepth();
	
	
	public  AlignmentPos flipped()
	{
		return this.clone(true);
	}
	
	
	
	/** These really only apply to SNPs */

	public Vector<ReadPos> getReadPositions()
	{
		return getReadPositions(false);
	}

	public Vector<ReadPos> getReadPositions(boolean fwOnly)
	{
		System.err.println("Base class AlignmentPos can not execute getReadPositions(fwOnly).  Use AlignmentPosSnps instead");
		(new Exception()).printStackTrace();
		return new Vector<ReadPos>();
	}
	
	public SymbolCounter getSnpCounter(boolean fwOnly)
	{
		System.err.println("Base class AlignmentPos can not execute getSnpCounter().  Use AlignmentPosSnps instead");
		(new Exception()).printStackTrace();
		return new SymbolCounter();
	}
	
	public SymbolCounterStratified getSnpCounterStratifiedByCycle(boolean fwOnly)
	{
		System.err.println("Base class AlignmentPos can not execute getSnpCounterStratified().  Use AlignmentPosSnps instead");
		(new Exception()).printStackTrace();
		return new SymbolCounterStratified();
	}

	public Vector<ReadPos> getReadPositions(StrandedFeature.Strand strand)
	{
		System.err.println("Base class AlignmentPos can not execute getReadPositions(strand).  Use AlignmentPosSnps instead");
		(new Exception()).printStackTrace();
		return new Vector<ReadPos>();
	}
	
	public SymbolCounter getSnpCounter(StrandedFeature.Strand strand)
	{
		System.err.println("Base class AlignmentPos can not execute getSnpCounter(strand).  Use AlignmentPosSnps instead");
		(new Exception()).printStackTrace();
		return new SymbolCounter();
	}
	
	public SymbolCounterStratified getSnpCounterStratifiedByCycle(StrandedFeature.Strand strand)
	{
		System.err.println("Base class AlignmentPos can not execute getSnpCounterStratified(strand).  Use AlignmentPosSnps instead");
		(new Exception()).printStackTrace();
		return new SymbolCounterStratified();
	}	
	
	
	
	
	
	public AlignmentPos clone()
	{
		return clone(false);
	}
	
	abstract public AlignmentPos clone(boolean flip_strand);


	
	

	/*****************
	 * 
	 * Setters
	 * 
	 */
	
	abstract public void removeRevStrandReads();
	abstract public void reset();
	

	
	/**
	 * @param coll
	 * @param fw
	 * @return
	 * 
	 * This function is destructive, so clone the collection first if you need
	 * a copy!
	 */
	static void	removeApsByStrand(Collection<? extends AlignmentPos> coll, boolean fw)
	{
		Iterator<? extends AlignmentPos> it = coll.iterator(); 
		while (it.hasNext())
		{
			AlignmentPos ap = it.next();
			if ((fw && (ap.getStrand() == StrandedFeature.POSITIVE)) ||
					(!fw && (ap.getStrand() == StrandedFeature.NEGATIVE)))
			{
				it.remove();
			}
		}
		
	}
	
	
	/** GFF **/
	
	public SimpleGFFRecord toGff(boolean flip)
	{
		AlignmentPos ap = this;
		if (flip)
		{
			ap = ap.flipped();
		}
		
		SimpleGFFRecord gff = new SimpleGFFRecord();
		gff.setStrand(ap.getStrand());
		gff.setFeature("exon");
		gff.setStart(ap.getPos());
		gff.setEnd(ap.getPos());
		gff.setSeqName(ap.getChr());
		gff.setSource("AlignmentPos");
		GFFUtils.add_gffrecord_map_entry(gff, "fw_depth", "" + ap.getDepth(true));
		GFFUtils.add_gffrecord_map_entry(gff, "rev_depth", "" + ap.getDepth(false));
		return gff;
	}
	
	
	public String gffLine()
	{
		return gffLine(false);
	}

	public String gffLine(boolean flip)
	{
		GFFRecord rec = this.toGff(flip);
		String out = GFFUtils.gffLine(rec);
		return out;
	}
	
	/*** Other output ***/
	@Override public String toString()
	{
		String out = null;

		out = this.chr + ":" + this.pos + this.strand.getToken() + "(" + Character.toUpperCase(this.getRefToken()) + ") depth= " +
		this.getDepth(true) + ", " + this.getDepth(false); 
		
		return out;
	}



}
