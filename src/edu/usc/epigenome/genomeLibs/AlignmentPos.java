package edu.usc.epigenome.genomeLibs;

import java.io.*;
import java.util.*;

import org.biojava.bio.program.gff.*;
import org.biojava.bio.seq.StrandedFeature;

public class AlignmentPos {

	/* Class vars */
	public static final double NO_COVERAGE = -1.0;

	/* Obj vars */
	public char f_ref;
	public String f_chr;
	public int f_pos;
	public int[] f_depth = null;



	/*****************
	 *  Constructors
	 */
	public AlignmentPos()
	{
	}
	
	public AlignmentPos(char ref, String chr, int pos)
	{
		f_ref = ref;
		f_chr = chr;
		f_pos = pos;
		f_depth = new int[] {0,0};
	}

	public AlignmentPos(char ref, String chr, int pos, int[] depth)
	{
		f_ref = ref;
		f_chr = chr;
		f_pos = pos;
		f_depth = depth;
	}
	
	
	/**
	 * 
	 * Getters
	 * 
	 */



	public char getRef(boolean reference_forward_strand)
	{
		char c = this.f_ref;
		if (!reference_forward_strand) c = MiscUtils.revCompNuc(c);
		return c;
	}
	
	
	public int getTotalDepth()
	{
		int[] d = getDepth();
		return d[0] + d[1];
	}
	
	
	public int[] getDepth()
	{
		return f_depth;
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
	
	public  AlignmentPos flipped ()
	{
		return this.clone(true);
	}
	
	public  AlignmentPos clone(boolean flip_strand)
	{
		char ref = (flip_strand) ? MiscUtils.revCompNuc(this.f_ref) : this.f_ref;
		AlignmentPos ap = new AlignmentPos(ref, this.f_chr, this.f_pos);

		if (f_depth != null)
		{
			ap.f_depth = new int[2];
			ap.f_depth[0] = this.f_depth[1];
			ap.f_depth[1] = this.f_depth[0];
		}
		
		return ap;
	}
	

	/*****************
	 * 
	 * Setters
	 * 
	 */
	
	
	public void add(char c, boolean read_forward_strand)
	throws Exception
	{
		add(read_forward_strand);
	}
	
	public void add(boolean read_forward_strand)
	throws Exception
	{
		int[] addition = new int[] {0,0};
		addition[ (read_forward_strand?0:1)]++;
		addDepth(addition);
	}

	public void addDepth(int[] depth)
	{
		f_depth[0] += depth[0];
		f_depth[1] += depth[1];
	}
	


	// Returns the actual nucleotide added
	public char addMaqPileupChar(char c)
	throws Exception
	{
		// Check for maq special characters
		char readstrand_c;
		boolean readstrand_forward_strand;
		
		switch(c)
		{
		case ',': readstrand_c = f_ref; readstrand_forward_strand = true; break;
		case '.': readstrand_c = f_ref; readstrand_forward_strand = false; break;
		default: readstrand_c = c; readstrand_forward_strand = Character.isUpperCase(c); break;
		}

		add(readstrand_forward_strand);
		return readstrand_c;
	}

	
	public void resetCounts()
	{
		this.f_depth = new int[] {0,0};
	}
	
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

		out = f_chr + ":" + f_pos + "(" + f_ref + ") depth= " +
		this.getDepth(true) + ", " + this.getDepth(false); 
		
		return out;
	}



}
