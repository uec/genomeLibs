package edu.usc.epigenome.genomeLibs;

import java.util.Collection;

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
	protected AlignmentPosOptions apOptions = null;
	


	/*****************
	 *  Constructors
	 */
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
		Symbol newRef = null;
		try
		{
			newRef = DNATools.complement(this.ref); 
		}
		catch (IllegalSymbolException e)
		{
			System.err.println(e);
			System.exit(0);
		}
		return newRef;
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
	
	/** GFF **/
	
	public SimpleGFFRecord toGff(String chrom, boolean ref_fw_strand, int pos)
	{
		SimpleGFFRecord gff = new SimpleGFFRecord();
		gff.setStrand((ref_fw_strand)? StrandedFeature.POSITIVE : StrandedFeature.NEGATIVE);
		gff.setFeature("exon");
		gff.setStart(pos);
		gff.setEnd(pos);
		gff.setSeqName(chrom);
		gff.setSource("AlignmentPos");
		GFFUtils.add_gffrecord_map_entry(gff, "fw_depth", "" + this.getDepth(true));
		GFFUtils.add_gffrecord_map_entry(gff, "rev_depth", "" + this.getDepth(false));
		return gff;
	}
	
	
	
	/*** Other output ***/
	@Override public String toString()
	{
		String out = null;

		out = this.chr + ":" + this.pos + "(" + Character.toUpperCase(this.getRefToken()) + ") depth= " +
		this.getDepth(true) + ", " + this.getDepth(false); 
		
		return out;
	}



}
