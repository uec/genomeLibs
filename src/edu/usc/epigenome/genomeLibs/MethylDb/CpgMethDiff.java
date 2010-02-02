package edu.usc.epigenome.genomeLibs.MethylDb;


/**
 * @author benb
 * A fake Cpg that just reports the methylation diff between two original Cpgs
 */
public class CpgMethDiff extends Cpg {
	
	protected Cpg cpga = null;
	protected Cpg cpgb = null;
	/**
	 * 
	 */
	public CpgMethDiff(Cpg inA, Cpg inB) {
		super(inA.chromPos, inA.negStrand);

		init(inA,inB);
	}
	
	
	protected void init(Cpg inA, Cpg inB)
	{
		cpga = inA;
		cpgb = inB;
	}


	/* (non-Javadoc)
	 * @see edu.usc.epigenome.genomeLibs.MethylDb.Cpg#fracMeth(boolean)
	 * 
	 * Returns a-b
	 */
	@Override
	public double fracMeth(boolean useNonconvFilt) 
	{

		double mA = cpga.fracMeth(useNonconvFilt);
		double mB = cpgb.fracMeth(useNonconvFilt);
		
		double out = mA-mB;
		
		return out;
	}
	
	

}
