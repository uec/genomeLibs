package edu.usc.epigenome.genomeLibs.MethylDb.CpgWalker;

import edu.usc.epigenome.genomeLibs.MethylDb.Cpg;

public class CpgWalkerAllpairsAutocorrByreadWcontext extends
		CpgWalkerAllpairsAutocorrByread {

	protected String fromContext = null;
	protected String toContext = null;
	
	public CpgWalkerAllpairsAutocorrByreadWcontext(
			CpgWalkerParams inWalkParams, boolean inSamestrandOnly,
			boolean inOppstrandOnly, boolean inSameRead, boolean inDifferentRead,
			String inFromContext, String inToContext) {
		super(inWalkParams, inSamestrandOnly, inOppstrandOnly, inSameRead,
				inDifferentRead);
		
		fromContext = inFromContext;
		toContext = inToContext;
	}

	@Override
	protected void recordPair(Cpg a, Cpg b) {
		
		
		String ac = a.context();
		String bc = b.context();
		
		if (ac.equals(this.fromContext) && bc.equals(this.toContext))
		{
			//System.err.printf("%s,%s context match (%s,%s)\n",a.chromPos, b.chromPos, ac, bc);
			super.recordPair(a, b);
		}
//		else
//		{
//			System.err.printf("%s,%s context MISMATCH (%s,%s) should be (%s,%s)\n",
//					a.chromPos, b.chromPos, ac, bc, fromContext, toContext);
//		}
	}
	
	
	

}
