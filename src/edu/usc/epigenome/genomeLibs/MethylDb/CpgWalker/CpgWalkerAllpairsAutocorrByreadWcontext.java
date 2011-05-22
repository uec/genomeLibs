package edu.usc.epigenome.genomeLibs.MethylDb.CpgWalker;

import edu.usc.epigenome.genomeLibs.MethylDb.Cpg;

public class CpgWalkerAllpairsAutocorrByreadWcontext extends
		CpgWalkerAllpairsAutocorrByread {

	protected String fromContext = null;
	protected String toContext = null;
	private int numPre = 0;
	private int numPost = 0;
	
	public CpgWalkerAllpairsAutocorrByreadWcontext(
			CpgWalkerParams inWalkParams, boolean inSamestrandOnly,
			boolean inOppstrandOnly, boolean inSameRead, boolean inDifferentRead,
			String inFromContext, String inToContext) {
		super(inWalkParams, inSamestrandOnly, inOppstrandOnly, inSameRead,
				inDifferentRead);
		
		fromContext = inFromContext;
		toContext = inToContext;
		
		if (fromContext.length() == 2)
		{
			numPre = 0;
			numPost = 1;
		}
		else if (fromContext.length() == 3)
		{
			numPre = 1;
			numPost = 1;
		}
	}

	@Override
	protected void recordPair(Cpg a, Cpg b) {
		
		
		String ac = a.context(numPre, numPost);
		String bc = b.context(numPre, numPost);
		
		if (ac.equals(this.fromContext) && bc.equals(this.toContext))
		{
			//System.err.printf("%s,%s context match (%s,%s)\n",a.chromPos, b.chromPos, ac, bc);
			super.recordPair(a, b);
		}
		else if (ac.equals(this.toContext) && bc.equals(this.fromContext))
		{
			super.recordPair(b, a);
		}

		
//		else
//		{
//			System.err.printf("%s,%s context MISMATCH (%s,%s) should be (%s,%s)\n",
//					a.chromPos, b.chromPos, ac, bc, fromContext, toContext);
//		}
	}
	
	
	

}
