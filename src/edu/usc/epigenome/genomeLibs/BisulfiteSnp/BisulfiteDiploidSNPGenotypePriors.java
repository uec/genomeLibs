package edu.usc.epigenome.genomeLibs.BisulfiteSnp;

import org.broadinstitute.sting.gatk.walkers.genotyper.DiploidSNPGenotypePriors;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.genotype.DiploidGenotype;

public class BisulfiteDiploidSNPGenotypePriors extends DiploidSNPGenotypePriors {

	protected static double Bisulfite_conversion_rate = 0.5;
	
	public BisulfiteDiploidSNPGenotypePriors() {
		// TODO Auto-generated constructor stub
	}

	public BisulfiteDiploidSNPGenotypePriors(byte ref, double heterozygosity,
			double probOfTriStateGenotype) {
		super(ref, heterozygosity, probOfTriStateGenotype);
		// TODO Auto-generated constructor stub
	}

	public BisulfiteDiploidSNPGenotypePriors(double[] log10Priors) {
		super(log10Priors);
		// TODO Auto-generated constructor stub
	}

	//@Override why can't be override??
	/**
     * Takes reference base, and three priors for hom-ref, het, hom-var, and fills in the priors vector
     * appropriately.
     *
     * Suppose A is the reference base, and we are given the probability of being hom-ref, het, and hom-var,
     * and that pTriSateGenotype is the true probability of observing reference A and a true genotype of B/C
     * then this sets the priors to:
     *
     * AA = hom-ref
     * AC = AG = AT = (het - pTriStateGenotype) / 3
     * CC = GG = TT = hom-var / 3
     * CG = CT = GT = pTriStateGenotype / 3
     *
     * So that we get:
     *
     * hom-ref + 3 * (het - pTriStateGenotype) / 3 + 3 * hom-var / 3 + 3 * pTriStateGenotype
     * hom-ref + het - pTriStateGenotype + hom-var + pTriStateGenotype
     * hom-ref + het + hom-var
     * = 1
     *
     * @param ref
     * @param heterozyosity
     * @param pRefError
     */
    public static double[] getReferencePolarizedPriors(byte ref, double heterozyosity, double pRefError ) {
        if ( ! MathUtils.isBounded(pRefError, 0.0, 0.01) ) {
            throw new RuntimeException(String.format("BUG: p Reference error is out of bounds (0.0 - 0.01) is allow range %f", pRefError));
        }

        double pTriStateGenotype = heterozyosity * pRefError;
//        if ( pTriStateGenotype >= heterozyosity ) {
//            throw new RuntimeException(String.format("p Tristate genotype %f is greater than the heterozygosity %f", pTriStateGenotype, heterozyosity));
//        }

        double pHomRef = heterozygosity2HomRefProbability(heterozyosity);
        double pHet    = heterozygosity2HetProbability(heterozyosity);
        double pHomVar = heterozygosity2HomVarProbability(heterozyosity);

        if (MathUtils.compareDoubles(pHomRef + pHet + pHomVar, 1.0) != 0) {
            throw new RuntimeException(String.format("BUG: Prior probabilities don't sum to one => %f, %f, %f", pHomRef, pHet, pHomVar));
        }

        double[] priors = new double[DiploidGenotype.values().length];

        for ( DiploidGenotype g : DiploidGenotype.values() ) {
            double POfG;

            final double nOnRefHets = 3;
            final double nOffRefHets = 3;
            final double nHomVars = 3;

            if ( g.isHomRef(ref) ){  
            	if(ref == 'C'|| ref == 'T'){
            		POfG = pHomRef*(1-Bisulfite_conversion_rate);
            	}
            	else{
            		POfG = pHomRef;
            	}
            }
            else if ( g.isHomVar(ref) ){ 
            	if((ref == 'C' && g.base1 == 'T')|| (ref == 'T' && g.base2 == 'C')){
            		POfG = pHomRef*Bisulfite_conversion_rate +  pHomVar/ nHomVars;
            	}
            	else{
            		POfG = pHomVar / nHomVars;
            	}
            }
            else if ( g.isHetRef(ref) ){ 
            	if(ref == 'C'|| ref == 'T'){
            		POfG = pHomRef*(1-Bisulfite_conversion_rate);
            	}
            	else{
            		POfG = (pHet - pTriStateGenotype ) / nOnRefHets;
            	}
            }
            else{ 
            	POfG = pTriStateGenotype / nOffRefHets; 
            }

            priors[g.ordinal()] = Math.log10(POfG);
        }

        return priors;
    }
    
}
