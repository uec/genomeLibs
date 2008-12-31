/**
 * 
 */
package edu.usc.epigenome.genomeLibs;


/**
 * @author benb
 * 
 *
 */
public class APHandlerWindowCounts extends GenomicWindCounter implements AlignmentPosStreamHandler {

	private static final long serialVersionUID = 8989385550284144139L;

	/**
	 * 
	 */
	public APHandlerWindowCounts(int inWindSize) {
		super(inWindSize);
	}

	/* (non-Javadoc)
	 * @see edu.usc.epigenome.genomeLibs.AlignmentPosStreamHandler#streamElement(java.util.LinkedList, edu.usc.epigenome.genomeLibs.AlignmentPos, java.util.LinkedList)
	 */
	public boolean streamElement(AlignmentPos[] pre, AlignmentPos currentAp, AlignmentPos[] post) 
	{
		this.increment(currentAp.getChr(), currentAp.getPos(), currentAp.getTotalDepth());
		return true;
	}

	/* (non-Javadoc)
	 * @see edu.usc.epigenome.genomeLibs.TreeMapCounter#excelOutput()
	 */
	@Override
	public String excelOutput() {
		// TODO Auto-generated method stub
		return super.excelOutput();
	}

	
	
	

}
