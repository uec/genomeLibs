/**
 * 
 */
package edu.usc.epigenome.genomeLibs;



import org.biojava.bio.seq.StrandedFeature;
import org.biojava.bio.symbol.Symbol;

/**
 * @author benb
 *
 */
public class ReadPosRich extends ReadPos {

	/* Obj vars */
	protected int cycle = UNKNOWN_CYCLE;
	protected int qual = UNKNOWN_QUAL;
	
	/** Constructors **/
	
	public ReadPosRich(Symbol inSym, boolean inForwardStrand, int inPos, int inQual)
	{
		this.sym = inSym;
		this.strand = (inForwardStrand) ? StrandedFeature.POSITIVE : StrandedFeature.NEGATIVE;
		this.qual = inQual;
		this.cycle = inPos;
	}

	public ReadPosRich(ReadPos inRp, int inPos, int inQual)
	{
		this.sym = inRp.getSym();
		this.strand = inRp.getStrand();
		this.qual = inQual;
		this.cycle = inPos;
	}
	
	
	/* Getters/Setters */

	/**
	 * @return the readPos
	 */
	public int getCycle() {
		return cycle;
	}


	/**
	 * @param readPos the readPos to set
	 */
	public void setReadPos(int pos) {
		this.cycle = pos;
	}


	/**
	 * @return the qual
	 */
	public int getQual() {
		return qual;
	}


	/**
	 * @param qual the qual to set
	 */
	public void setQual(int qual) {
		this.qual = qual;
	}
	
	
	
	
	
	
	
	
}
