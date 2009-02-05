package edu.usc.epigenome.genomeLibs.AlignmentPos;

import org.biojava.bio.symbol.Symbol;


public class AlignmentPosNull extends AlignmentPos {

	public AlignmentPosNull() {
		super('N',AlignmentPos.NULL_CHROM, -1, new AlignmentPosOptions());
	}

	@Override
	public AlignmentPos clone(boolean flip_strand) {
		return new AlignmentPosNull();
	}

	@Override
	public int[] getDepth() {
		return new int[] {0,0};
	}

	@Override
	public void reset() {
	}

	/* (non-Javadoc)
	 * @see edu.usc.epigenome.genomeLibs.AlignmentPos#setApOptions(edu.usc.epigenome.genomeLibs.AlignmentPosOptions)
	 */
	@Override
	public void setApOptions(AlignmentPosOptions inApOptions) {
		System.err.println("Can not set fields for an AlignmentPosNull obect, quitting");
		new Exception().printStackTrace();
		System.exit(0);
	}
	

	/* (non-Javadoc)
	 * @see edu.usc.epigenome.genomeLibs.AlignmentPos#setChr(java.lang.String)
	 */
	@Override
	public void setChr(String chr) {
		System.err.println("Can not set fields for an AlignmentPosNull obect, quitting");
		new Exception().printStackTrace();
		System.exit(0);
	}

	/* (non-Javadoc)
	 * @see edu.usc.epigenome.genomeLibs.AlignmentPos#setPos(int)
	 */
	@Override
	public void setPos(int pos) {
		System.err.println("Can not set fields for an AlignmentPosNull obect, quitting");
		new Exception().printStackTrace();
		System.exit(0);
	}

	/* (non-Javadoc)
	 * @see edu.usc.epigenome.genomeLibs.AlignmentPos#setRef(org.biojava.bio.symbol.Symbol)
	 */
	@Override
	public void setRef(Symbol inRef) {
		System.err.println("Can not set fields for an AlignmentPosNull obect, quitting");
		new Exception().printStackTrace();
		System.exit(0);
	}

	/* (non-Javadoc)
	 * @see edu.usc.epigenome.genomeLibs.AlignmentPos#removeRevStrandReads()
	 */
	@Override
	public void removeRevStrandReads() {}


}
