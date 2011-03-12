package org.broadinstitute.sting.gatk.bisulfitegenotyper;

import java.util.Arrays;

import org.broad.tribble.util.variantcontext.VariantContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.utils.helpers.DbSNPHelper;
import org.broadinstitute.sting.gatk.walkers.genotyper.GenotypePriors;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.genotype.DiploidGenotype;

public class BisulfiteDiploidSNPGenotypePriors implements GenotypePriors {
	// --------------------------------------------------------------------------------------------------------------
    //
    // Constants and static information
    //
    // --------------------------------------------------------------------------------------------------------------
    public static final double HUMAN_HETEROZYGOSITY = 1e-3;
    public static final double CEU_HETEROZYGOSITY = 1e-3;
    public static final double YRI_HETEROZYGOSITY = 1.0 / 850;
    public static final double DBSNP_HETEROZYGOSITY = 0.5;

    
    /**
     * Default value of the prob of seeing a reference error.  Used to calculation the
     * chance of seeing a true B/C het when the reference is A, which we assume is the product
     * of the ref error rate and the het. Default value is Q60
     */
    public static final double PROB_OF_REFERENCE_ERROR = 1e-6;  // the reference is
    
    protected static double Bisulfite_conversion_rate = 0.95;
    
 
    private final static double[] flatPriors = new double[DiploidGenotype.values().length];

    // --------------------------------------------------------------------------------------------------------------
    //
    // Diploid priors
    //
    // --------------------------------------------------------------------------------------------------------------
    private double[] priors = null;

    // todo -- fix me when this issue is resolved
    public static final boolean requirePriorSumToOne = false;

    /**
     * Create a new DiploidGenotypePriors object with flat priors for each diploid genotype
     */
    public BisulfiteDiploidSNPGenotypePriors() {
        priors = flatPriors.clone();
    }

    /**
     * Create a new GenotypeLikelihoods object with priors for a diploid with heterozygosity and reference
     * base ref
     *
     * @param ref
     * @param heterozygosity
     * @param probOfTriStateGenotype The prob of seeing a true B/C het when the reference is A
     */
    public BisulfiteDiploidSNPGenotypePriors(byte ref, double heterozygosity, double probOfTriStateGenotype) {
        priors = getReferencePolarizedPriors(ref, heterozygosity, probOfTriStateGenotype);
    }
    
    public BisulfiteDiploidSNPGenotypePriors(byte ref, double heterozygosity, double probOfTriStateGenotype, double bisulfiteConversionRate) {
    	Bisulfite_conversion_rate = bisulfiteConversionRate;
    	priors = getReferencePolarizedPriors(ref, heterozygosity, probOfTriStateGenotype);
        
    }

    /**
     * Create a new Genotypelike Likelhoods's object with priors (in log10 space) for each of the DiploteGenotypes
     *
     * @param log10Priors
     */
    public BisulfiteDiploidSNPGenotypePriors(double[] log10Priors) {
        priors = log10Priors.clone();
    }

    /**
     * Returns an array of priors for each genotype, indexed by DiploidGenotype.ordinal values().
     *
     * @return log10 prior as a double array
     */
    public double[] getPriors() {
        return priors;
    }
    
    public void setPriors(RefMetaDataTracker tracker, ReferenceContext ref, double heterozygosity, double probOfTriStateGenotype, double bisulfiteConversionRate) {
    	Bisulfite_conversion_rate = bisulfiteConversionRate;
    	String rsID = DbSNPHelper.rsIDOfFirstRealSNP(tracker.getReferenceMetaData(DbSNPHelper.STANDARD_DBSNP_TRACK_NAME));
    	if ( rsID != null ){
    		priors = getReferencePolarizedPriors(ref.getBase(), DBSNP_HETEROZYGOSITY, probOfTriStateGenotype);
    	}
    	else{
    		priors = getReferencePolarizedPriors(ref.getBase(), heterozygosity, probOfTriStateGenotype);
    	}
    	
    }

    /**
     * Returns the prior associated with DiploidGenotype g
     * @param g
     * @return log10 prior as a double
     */
    public double getPrior(DiploidGenotype g) {
        return getPriors()[g.ordinal()];
    }

    public double getHeterozygosity() { return HUMAN_HETEROZYGOSITY; }

    public boolean validate(boolean throwException) {
        try {
            if ( requirePriorSumToOne && MathUtils.compareDoubles(MathUtils.sumLog10(priors), 1.0) != 0 ) {
                throw new IllegalStateException(String.format("Priors don't sum to 1: sum=%f %s", MathUtils.sumLog10(priors), Arrays.toString(priors)));
            }

            for ( DiploidGenotype g : DiploidGenotype.values() ) {
                int i = g.ordinal();
                if ( ! MathUtils.wellFormedDouble(priors[i]) || ! MathUtils.isNegativeOrZero(priors[i]) ) {
                    String bad = String.format("Prior %f is badly formed %b", priors[i], MathUtils.isNegativeOrZero(priors[i]));
                    throw new IllegalStateException(String.format("At %s: %s", g.toString(), bad));
                }
            }
        } catch ( IllegalStateException e ) {
            if ( throwException )
                throw new RuntimeException(e);
            else
                return false;
        }

        return true;
    }

    // --------------------------------------------------------------------------------------------------------------
    //
    // Static functionality
    //
    // --------------------------------------------------------------------------------------------------------------

    /**
     * Returns homozygous-reference, heterozygous, and homozygous-non-ref probabilities given a heterozygosity
     * value, as elements 0, 1, and 2 of a double[], respectively
     *
     * @param h the heterozygosity [probability of a base being heterozygous]
     */
    @Deprecated
    public static double[] heterozygosity2DiploidProbabilities(double h) {
        double[] pdbls = new double[3];

        pdbls[0] = heterozygosity2HomRefProbability(h);
        pdbls[1] = heterozygosity2HetProbability(h);
        pdbls[2] = heterozygosity2HomVarProbability(h);
        return pdbls;
    }

    /**
     *
     * @param h
     * @return
     */
    public static double heterozygosity2HomRefProbability(double h) {
        if (MathUtils.isNegative(h)) {
            throw new RuntimeException(String.format("Heterozygous value is bad %f", h));
        }

        double v = 1.0 - (3.0 * h / 2.0);
        if (MathUtils.isNegative(v)) {
            throw new RuntimeException(String.format("Heterozygous value is bad %f", h));
        }

        return v;
    }

    public static double heterozygosity2HetProbability(double h) {
        if (MathUtils.isNegative(h)) {
            throw new RuntimeException(String.format("Heterozygous value is bad %f", h));
        }

        return h;
    }

    public static double heterozygosity2HomVarProbability(double h) {
        if (MathUtils.isNegative(h)) {
            throw new RuntimeException(String.format("Heterozygous value is bad %f", h));
        }

        return h / 2.0;
    }


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
     * when ref=A
     * AA = hom-ref
     * AG = (het - pTriStateGenotype) / 3
     * AC = ((het - pTriStateGenotype) / 3) * (1 - Bisulfite_conversion_rate)
     * AT = ((het - pTriStateGenotype) / 3) * (1 + Bisulfite_conversion_rate)
     * GG = hom-var / 3
     * CC = (hom-var / 3) * (1 - Bisulfite_conversion_rate)
     * TT = hom-var / 3 + (pTriStateGenotype / 3) * Bisulfite_conversion_rate
     * CG = (pTriStateGenotype / 3) * (1 - Bisulfite_conversion_rate)
     * GT = pTriStateGenotype / 3
     * CT = (pTriStateGenotype / 3) * (1 - Bisulfite_conversion_rate) + (hom-var / 3) * Bisulfite_conversion_rate
     * 
     * when ref=C
     * CC = hom-ref * (1 - Bisulfite_conversion_rate)
     * CA = CG = ((het - pTriStateGenotype) / 3) * (1 - Bisulfite_conversion_rate)
     * CT = ((het - pTriStateGenotype) / 3) * (1 - Bisulfite_conversion_rate) + hom-ref * Bisulfite_conversion_rate
     * AA = GG = hom-var / 3
     * TT = hom-var / 3 + ((het - pTriStateGenotype) / 3) * Bisulfite_conversion_rate
     * AG = pTriStateGenotype / 3
     * AT = GT = pTriStateGenotype / 3 + ((het - pTriStateGenotype) / 3) * Bisulfite_conversion_rate
     * 
     * when ref=T
     * TT = hom-ref + ((het - pTriStateGenotype) / 3) * Bisulfite_conversion_rate
     * TA = TG = (het - pTriStateGenotype) / 3 + (pTriStateGenotype / 3) * Bisulfite_conversion_rate 
     * CT = ((het - pTriStateGenotype) / 3) * (1 - Bisulfite_conversion_rate) + (hom-var / 3) * Bisulfite_conversion_rate
     * AA = GG = hom-var / 3
     * CC = (hom-var / 3) * (1 - Bisulfite_conversion_rate)
     * AG = pTriStateGenotype / 3
     * AC = GC = (pTriStateGenotype / 3) * (1 - Bisulfite_conversion_rate)
     * 
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
            	if(ref == BaseUtils.C){
            		POfG = pHomRef*(1-Bisulfite_conversion_rate);
            	}
            	else if(ref == BaseUtils.T){
            		POfG = pHomRef + (pHet - pTriStateGenotype) * Bisulfite_conversion_rate/3.0;
            	}
            	else{
            		POfG = pHomRef;
            	}
            }
            else if ( g.isHomVar(ref) ){ 
            	if(ref == BaseUtils.C && g.base1 == BaseUtils.T){
            		POfG = pHomVar/3.0 + (pHet - pTriStateGenotype) * Bisulfite_conversion_rate/3.0;
            	}
            	else if((ref == BaseUtils.A || ref == BaseUtils.G || ref == BaseUtils.T) && g.base1 == BaseUtils.C){
            		POfG = (pHomVar/3.0) * (1.0 - Bisulfite_conversion_rate);
            	}
            	else if((ref == BaseUtils.A || ref == BaseUtils.G) && g.base1 == BaseUtils.T){
            		 POfG = pHomVar/3.0 + (pTriStateGenotype / 3.0) * Bisulfite_conversion_rate;
            	}
            	else{
            		POfG = pHomVar / nHomVars;
            	}
            }
            else if ( g.isHetRef(ref) ){ 
            	if(ref == BaseUtils.C){
            	     if((g.base1 == BaseUtils.C && g.base2 == BaseUtils.T)|| (g.base1 == BaseUtils.T && g.base2 == BaseUtils.C)){
            	    	 POfG = ((pHet - pTriStateGenotype) / 3.0) * (1 - Bisulfite_conversion_rate) + pHomRef * Bisulfite_conversion_rate;
            	     }
            	     else{
            	    	 POfG = ((pHet - pTriStateGenotype) / 3.0) * (1 - Bisulfite_conversion_rate);
            	     }
            	}
            	else if(ref == BaseUtils.T){
            		if((g.base1 == BaseUtils.C && g.base2 == BaseUtils.T)|| (g.base1 == BaseUtils.T && g.base2 == BaseUtils.C)){
           	    	 	POfG = ((pHet - pTriStateGenotype) / 3.0) * (1 - Bisulfite_conversion_rate) + (pHomVar/3.0) * Bisulfite_conversion_rate;
           	     	}
           	     	else{
           	     		POfG = (pHet - pTriStateGenotype) / 3.0 + (pTriStateGenotype / 3) * Bisulfite_conversion_rate;
           	     	}
            	}
            	else{
            	     if((g.base1 == ref && g.base2 == BaseUtils.C)|| (g.base1 == BaseUtils.C && g.base2 == ref)){
            	    	 POfG = ((pHet - pTriStateGenotype ) / 3.0) * (1 - Bisulfite_conversion_rate);
            	     }
            	     else if((g.base1 == ref && g.base2 == BaseUtils.T)|| (g.base1 == BaseUtils.T && g.base2 == ref)){
            	    	 POfG = ((pHet - pTriStateGenotype ) / 3.0) * (1 + Bisulfite_conversion_rate);
            	     }
            	     else{
            	    	 POfG = (pHet - pTriStateGenotype ) / 3.0;
            	     }
            	}
            }
            else{ 
            	if(ref == BaseUtils.C){
            	     if((g.base1 == BaseUtils.A && g.base2 == BaseUtils.G)|| (g.base1 == BaseUtils.G && g.base2 == BaseUtils.A)){
            	    	 POfG = pTriStateGenotype / 3.0;
            	     }
            	     else{
            	    	 POfG = pTriStateGenotype / 3.0 + ((pHet - pTriStateGenotype) / 3.0) * Bisulfite_conversion_rate;
            	     }
            	}
            	else if(ref == BaseUtils.T){
            	     if((g.base1 == BaseUtils.A && g.base2 == BaseUtils.G)|| (g.base1 == BaseUtils.G && g.base2 == BaseUtils.A)){
            	    	 POfG = pTriStateGenotype / 3.0;
            	     }
            	     else{
            	    	 POfG = (pTriStateGenotype / 3.0) * (1 - Bisulfite_conversion_rate);
            	     } 
            	}
            	else{
            	     if((g.base1 == BaseUtils.C && g.base2 == BaseUtils.G)|| (g.base1 == BaseUtils.G && g.base2 == BaseUtils.C)){
            	    	 POfG =  (pTriStateGenotype / 3.0) * (1 - Bisulfite_conversion_rate);
            	     }
            	     else if((g.base1 == BaseUtils.G && g.base2 == BaseUtils.T)|| (g.base1 == BaseUtils.T && g.base2 == BaseUtils.G)){
            	    	 POfG = pTriStateGenotype / 3.0;
            	     }
            	     else{
            	    	 POfG =  (pTriStateGenotype / 3.0) * (1 - Bisulfite_conversion_rate) + (pHomVar/3.0) * Bisulfite_conversion_rate;
            	     }
            	}
            	 
            }

            priors[g.ordinal()] = Math.log10(POfG);
        }

        return priors;
    }

    static {
        for ( DiploidGenotype g : DiploidGenotype.values() ) {
            flatPriors[g.ordinal()] = Math.log10(1.0 / DiploidGenotype.values().length);
        }
    }

}
