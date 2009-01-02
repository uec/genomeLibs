package edu.usc.epigenome.genomeLibs;

import java.util.*;

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
	 * @param strand the strand to set
	 */
	public void setStrand(StrandedFeature.Strand strand) {
		this.strand = strand;
	}

	static String getRefTokens(Collection<AlignmentPos> aps)
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
		
	static String getRefTokens(AlignmentPos[] aps)
	{
		StringBuffer buf = new StringBuffer(aps.length);
		for (int i = 0; i < aps.length; i++)
		{
			buf.append(aps[i].getRefToken());
		}
		
		return buf.toString();
	}

	static String getPosString(Collection<AlignmentPos> aps)
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
