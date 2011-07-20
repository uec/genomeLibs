package org.broadinstitute.sting.gatk.uscec.bisulfitesnpmodel;

import java.util.Arrays;

import org.broad.tribble.Feature;
import org.broad.tribble.dbsnp.DbSNPFeature;
import org.broad.tribble.util.variantcontext.VariantContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.utils.GATKFeature;
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
    public static double DBSNP_VALIDATE_HETEROZYGOSITY = 0.1;
    public static double DBSNP_NOVAL_HETEROZYGOSITY = 0.02;

    
    /**
     * Default value of the prob of seeing a reference error.  Used to calculation the
     * chance of seeing a true B/C het when the reference is A, which we assume is the product
     * of the ref error rate and the het. Default value is Q60
     */
    public static final double PROB_OF_REFERENCE_ERROR = 1e-6;  // the reference is
    
    protected static double BISULFITE_CONVERSION_RATE = 0; // it is the wrong concept here. it is means Cytosine methylation rate here, not means the bisulfite conversion rate!!!
    protected static double CYTOSINE_METHYLATION_RATE = 0;
    
 
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
       // priors = getReferencePolarizedPriors(ref, heterozygosity, probOfTriStateGenotype);
    	priors = flatPriors.clone();
    }
    
    public BisulfiteDiploidSNPGenotypePriors(byte ref, double heterozygosity, double probOfTriStateGenotype, double bisulfiteConversionRate) {
    //	BISULFITE_CONVERSION_RATE = bisulfiteConversionRate;
    //	priors = getReferencePolarizedPriors(ref, heterozygosity, probOfTriStateGenotype);
    	priors = flatPriors.clone();
        
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
    
    public void setPriors(RefMetaDataTracker tracker, ReferenceContext ref, double heterozygosity, double probOfTriStateGenotype, double bisulfiteConversionRate, double cytosineMethyRate, double novelDbsnpHet, double validateDbsnpHet, CytosineTypeStatus cts) {
    	BISULFITE_CONVERSION_RATE = bisulfiteConversionRate;
    	CYTOSINE_METHYLATION_RATE = cytosineMethyRate;
    	DBSNP_NOVAL_HETEROZYGOSITY = novelDbsnpHet;
    	DBSNP_VALIDATE_HETEROZYGOSITY = validateDbsnpHet;
        //System.err.println(refWindow.length);
    	
        byte refBase = ref.getBase();

        
        
        DbSNPFeature d = DbSNPHelper.getFirstRealSNP(tracker.getReferenceMetaData(DbSNPHelper.STANDARD_DBSNP_TRACK_NAME));
        String rsID = DbSNPHelper.rsIDOfFirstRealSNP(tracker.getReferenceMetaData(DbSNPHelper.STANDARD_DBSNP_TRACK_NAME));
        if(rsID != null){
        	//System.err.println("is dbsnp");
        	if(d.getValidationStatus().equalsIgnoreCase("unknown"))
        		priors = getReferencePolarizedPriorsBasedOnMethyStatus(refBase, DBSNP_NOVAL_HETEROZYGOSITY, probOfTriStateGenotype);
        	else
        		priors = getReferencePolarizedPriorsBasedOnMethyStatus(refBase, DBSNP_VALIDATE_HETEROZYGOSITY, probOfTriStateGenotype);
        }
        else{
        	priors = getReferencePolarizedPriorsBasedOnMethyStatus(refBase, heterozygosity, probOfTriStateGenotype);
        }
    	//String rsID = DbSNPHelper.rsIDOfFirstRealSNP(tracker.getReferenceMetaData(DbSNPHelper.STANDARD_DBSNP_TRACK_NAME));
    	//if ( rsID != null ){
    		
    		//priors = getReferencePolarizedPriors(ref.getBase(), DBSNP_HETEROZYGOSITY, probOfTriStateGenotype);
    	//}
    	//else{
    		
    		//priors = getReferencePolarizedPriors(ref.getBase(), heterozygosity, probOfTriStateGenotype);
    	//}
    	
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

        double v = 1.0 - (3 * h / 2.0);
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
    
    public static boolean isTransition(byte base, byte ref) {
        boolean transition;
    	switch(ref){
        	case BaseUtils.A: transition = base == BaseUtils.G ? true : false; break;
        	case BaseUtils.G: transition = base == BaseUtils.A ? true : false; break;
        	case BaseUtils.C: transition = base == BaseUtils.T ? true : false; break;
        	case BaseUtils.T: transition = base == BaseUtils.C ? true : false; break;
        	default: throw new RuntimeException(String.format("ref base is bad "));
        }
    	return transition;
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
    /*
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

            final double transitionRate = 2.0/3.0;
            final double transversionRate = 1.0/6.0;
            //final double nHomVars = 3;

            if ( g.isHomRef(ref) ){  
            	if(ref == BaseUtils.C){
            		POfG = pHomRef*(1-BISULFITE_CONVERSION_RATE);
            	}
            	else if(ref == BaseUtils.T){

            			POfG = pHomRef + (pHet - pTriStateGenotype) * BISULFITE_CONVERSION_RATE * transitionRate;
            		
            		
            	}
            	else{
            		POfG = pHomRef;
            	}
            }
            else if ( g.isHomVar(ref) ){ 
            	if(ref == BaseUtils.C && g.base1 == BaseUtils.T){
            		POfG = pHomVar * transitionRate + (pHet - pTriStateGenotype) * BISULFITE_CONVERSION_RATE * transitionRate;
            	}
            	else if((ref == BaseUtils.A || ref == BaseUtils.G || ref == BaseUtils.T) && g.base1 == BaseUtils.C){
            		
            		if(isTransition(g.base1,ref)){
            			POfG = (pHomVar * transitionRate) * (1.0 - BISULFITE_CONVERSION_RATE);
            		}
            		else{
            			POfG = (pHomVar * transversionRate) * (1.0 - BISULFITE_CONVERSION_RATE);
            		}
            		
            	}
            	else if((ref == BaseUtils.A || ref == BaseUtils.G) && g.base1 == BaseUtils.T){
            		 POfG = pHomVar * transversionRate + (pTriStateGenotype * transversionRate) * BISULFITE_CONVERSION_RATE;
            	}
            	else{
            		if(isTransition(g.base1,ref)){
            			POfG = pHomVar * transitionRate;
            		}
            		else{
            			POfG = pHomVar * transversionRate;
            		}
            		
            	}
            }
            else if ( g.isHetRef(ref) ){ 
            	if(ref == BaseUtils.C){
            	     if((g.base1 == BaseUtils.C && g.base2 == BaseUtils.T)|| (g.base1 == BaseUtils.T && g.base2 == BaseUtils.C)){
            	    	 POfG = ((pHet - pTriStateGenotype) * transitionRate) * (1 - BISULFITE_CONVERSION_RATE) + pHomRef * BISULFITE_CONVERSION_RATE;
            	     }
            	     else{
            	    	 POfG = ((pHet - pTriStateGenotype) * transversionRate) * (1 - BISULFITE_CONVERSION_RATE);
            	     }
            	}
            	else if(ref == BaseUtils.T){
            		if((g.base1 == BaseUtils.C && g.base2 == BaseUtils.T)|| (g.base1 == BaseUtils.T && g.base2 == BaseUtils.C)){
           	    	 	POfG = ((pHet - pTriStateGenotype) * transitionRate) * (1 - BISULFITE_CONVERSION_RATE) + (pHomVar * transitionRate) * BISULFITE_CONVERSION_RATE;
           	     	}
           	     	else{
           	     		POfG = (pHet - pTriStateGenotype) * transversionRate + (pTriStateGenotype * transversionRate) * BISULFITE_CONVERSION_RATE;
           	     	}
            	}
            	else{
            	     if((g.base1 == ref && g.base2 == BaseUtils.C)|| (g.base1 == BaseUtils.C && g.base2 == ref)){
            	    	 POfG = ((pHet - pTriStateGenotype ) * transversionRate) * (1 - BISULFITE_CONVERSION_RATE);
            	     }
            	     else if((g.base1 == ref && g.base2 == BaseUtils.T)|| (g.base1 == BaseUtils.T && g.base2 == ref)){
            	    	 POfG = ((pHet - pTriStateGenotype ) * transversionRate) * (1 + BISULFITE_CONVERSION_RATE);
            	     }
            	     else{
            	    	 POfG = (pHet - pTriStateGenotype ) * transitionRate;
            	     }
            	}
            }
            else{ 
            	if(ref == BaseUtils.C){
            	     if((g.base1 == BaseUtils.A && g.base2 == BaseUtils.G)|| (g.base1 == BaseUtils.G && g.base2 == BaseUtils.A)){
            	    	 POfG = pTriStateGenotype * transversionRate;
            	     }
            	     else{
            	    	 POfG = pTriStateGenotype * transitionRate + ((pHet - pTriStateGenotype) * transitionRate) * BISULFITE_CONVERSION_RATE;
            	     }
            	}
            	else if(ref == BaseUtils.T){
            	     if((g.base1 == BaseUtils.A && g.base2 == BaseUtils.G)|| (g.base1 == BaseUtils.G && g.base2 == BaseUtils.A)){
            	    	 POfG = pTriStateGenotype * transversionRate;
            	     }
            	     else{
            	    	 POfG = (pTriStateGenotype * transitionRate) * (1 - BISULFITE_CONVERSION_RATE);
            	     } 
            	}
            	else{
            	     if((g.base1 == BaseUtils.C && g.base2 == BaseUtils.G)|| (g.base1 == BaseUtils.G && g.base2 == BaseUtils.C)){
            	    	 if(isTransition(g.base1,ref) || isTransition(g.base2,ref)){
            	    		 POfG =  (pTriStateGenotype * transitionRate) * (1 - BISULFITE_CONVERSION_RATE);
            	    	 }
            	    	 else{
            	    		 POfG =  (pTriStateGenotype * transversionRate) * (1 - BISULFITE_CONVERSION_RATE);
            	    	 }
            	    	 
            	     }
            	     else if((g.base1 == BaseUtils.G && g.base2 == BaseUtils.T)|| (g.base1 == BaseUtils.T && g.base2 == BaseUtils.G)){
            	    	 if(isTransition(g.base1,ref) || isTransition(g.base2,ref)){
            	    		 POfG = pTriStateGenotype * transitionRate;
            	    	 }
            	    	 else{
            	    		 POfG = pTriStateGenotype * transversionRate;
            	    	 }
            	    	 
            	     }
            	     else{
            	    	 if(isTransition(g.base1,ref) || isTransition(g.base2,ref)){
            	    		 POfG =  (pTriStateGenotype * transitionRate) * (1 - BISULFITE_CONVERSION_RATE) + (pHomVar * transitionRate) * BISULFITE_CONVERSION_RATE;
            	    	 }
            	    	 else{
            	    		 POfG =  (pTriStateGenotype * transversionRate) * (1 - BISULFITE_CONVERSION_RATE) + (pHomVar * transversionRate) * BISULFITE_CONVERSION_RATE;
            	    	 }
            	    	 
            	     }
            	}
            	 
            }

            priors[g.ordinal()] = Math.log10(POfG);
        }

        return priors;
    }
    */
    
    
    /**
     * Takes reference base, and three priors for hom-ref, het, hom-var, and fills in the priors vector
     * appropriately.
     *
     * Suppose A is the reference base, and we are given the probability of being hom-ref, het, and hom-var,
     * and that pTriSateGenotype is the true probability of observing reference A and a true genotype of B/C
     * then this sets the priors to:
     *
     * AA = hom-ref
     * AC = AT = (het - pTriStateGenotype) * transversionRate
     * AG = (het - pTriStateGenotype) * transitionRate
     * CC = TT = hom-var * transversionRate
     * GG = hom-var * transitionRate
     * CG = GT = pTriStateGenotype * transitionRate
     * CT = pTriStateGenotype * transversionRate
     * 
     * when ref=A
     * AA = hom-ref
     * AG = (het - pTriStateGenotype) * transitionRate
     * AC = ((het - pTriStateGenotype) * transversionRate) * ((1 - Cytosine_methylation_rate) * (1 - Bisulfite_conversion_rate) + Cytosine_methylation_rate)
     * AT = ((het - pTriStateGenotype) * transversionRate) * (1 + (1 - Cytosine_methylation_rate) * Bisulfite_conversion_rate)
     * GG = hom-var * transitionRate
     * CC = (hom-var * transversionRate) * ((1 - Cytosine_methylation_rate) * (1 - Bisulfite_conversion_rate) + Cytosine_methylation_rate)
     * TT = hom-var * transversionRate + (pTriStateGenotype * transversionRate) * (1 - Cytosine_methylation_rate) * Bisulfite_conversion_rate
     * CG = (pTriStateGenotype * transitionRate) * ((1 - Cytosine_methylation_rate) * (1 - Bisulfite_conversion_rate) + Cytosine_methylation_rate)
     * GT = pTriStateGenotype * transitionRate
     * CT = (pTriStateGenotype * transversionRate) * ((1 - Cytosine_methylation_rate) * (1 - Bisulfite_conversion_rate) + Cytosine_methylation_rate) + (hom-var * transversionRate) * (1 - Cytosine_methylation_rate) * Bisulfite_conversion_rate
     * 
     * when ref=C
     * CC = hom-ref * ((1 - Cytosine_methylation_rate) * (1 - Bisulfite_conversion_rate) + Cytosine_methylation_rate)
     * CA = CG = ((het - pTriStateGenotype) * transversionRate) * ((1 - Cytosine_methylation_rate) * (1 - Bisulfite_conversion_rate) + Cytosine_methylation_rate)
     * CT = ((het - pTriStateGenotype) * transitionRate) * ((1 - Cytosine_methylation_rate) * (1 - Bisulfite_conversion_rate) + Cytosine_methylation_rate) + hom-ref * (1 - Cytosine_methylation_rate) * Bisulfite_conversion_rate
     * AA = GG = hom-var * transversionRate
     * TT = hom-var * transitionRate + ((het - pTriStateGenotype) * transitionRate) * (1 - Cytosine_methylation_rate) * Bisulfite_conversion_rate
     * AG = pTriStateGenotype * transversionRate
     * AT = GT = pTriStateGenotype * transitionRate + ((het - pTriStateGenotype) * transitionRate) * (1 - Cytosine_methylation_rate) * Bisulfite_conversion_rate
     * 
     * when ref=T
     * TT = hom-ref + ((het - pTriStateGenotype) * transitionRate) * (1 - Cytosine_methylation_rate) * Bisulfite_conversion_rate
     * TA = TG = (het - pTriStateGenotype) * transversionRate + (pTriStateGenotype * transversionRate) * (1 - Cytosine_methylation_rate) * Bisulfite_conversion_rate 
     * CT = ((het - pTriStateGenotype) * transitionRate) * ((1 - Cytosine_methylation_rate) * (1 - Bisulfite_conversion_rate) + Cytosine_methylation_rate) + (hom-var * transitionRate) * (1 - Cytosine_methylation_rate) * Bisulfite_conversion_rate
     * AA = GG = hom-var * transversionRate
     * CC = (hom-var * transitionRate) * ((1 - Cytosine_methylation_rate) * (1 - Bisulfite_conversion_rate) + Cytosine_methylation_rate)
     * AG = pTriStateGenotype * transversionRate
     * AC = GC = (pTriStateGenotype * transitionRate) * ((1 - Cytosine_methylation_rate) * (1 - Bisulfite_conversion_rate) + Cytosine_methylation_rate)
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
    
    public static double[] getReferencePolarizedPriorsBasedOnMethyStatus(byte ref, double heterozyosity, double pRefError) {
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
            double POfG = 1;

            final double transitionRate = 2.0/3.0;
            final double transversionRate = 1.0/6.0;
            //final double nHomVars = 3;
            
            if(ref == BaseUtils.A){
            	if(g.base1 == BaseUtils.A){
            		if(g.base2 == BaseUtils.A){
            			POfG = pHomRef;
                    }
                    else if(g.base2 == BaseUtils.G){
                    	POfG = (pHet - pTriStateGenotype) * transitionRate;
                    }
                    else if(g.base2 == BaseUtils.C){
                    	POfG = ((pHet - pTriStateGenotype) * transversionRate) * ((1 - CYTOSINE_METHYLATION_RATE) * (1 - BISULFITE_CONVERSION_RATE) + CYTOSINE_METHYLATION_RATE);
                    }
                    else if(g.base2 == BaseUtils.T){
                    	POfG = ((pHet - pTriStateGenotype) * transversionRate) * (1 + (1 - CYTOSINE_METHYLATION_RATE) * BISULFITE_CONVERSION_RATE);
                    }
                    else{
                    	
                    }
                }
                else if(g.base1 == BaseUtils.G){
                	if(g.base2 == BaseUtils.A){
                		POfG = (pHet - pTriStateGenotype) * transitionRate;
                    }
                    else if(g.base2 == BaseUtils.G){
                    	POfG = pHomVar * transitionRate;
                    }
                    else if(g.base2 == BaseUtils.C){
                    	POfG = (pTriStateGenotype * transitionRate) * ((1 - CYTOSINE_METHYLATION_RATE) * (1 - BISULFITE_CONVERSION_RATE) + CYTOSINE_METHYLATION_RATE);
                    }
                    else if(g.base2 == BaseUtils.T){
                    	POfG = pTriStateGenotype * transitionRate;
                    }
                    else{
                    	
                    }
                }
                else if(g.base1 == BaseUtils.C){
                	if(g.base2 == BaseUtils.A){
                		POfG = ((pHet - pTriStateGenotype) * transversionRate) * ((1 - CYTOSINE_METHYLATION_RATE) * (1 - BISULFITE_CONVERSION_RATE) + CYTOSINE_METHYLATION_RATE);
                    }
                    else if(g.base2 == BaseUtils.G){
                    	POfG = (pTriStateGenotype * transitionRate) * ((1 - CYTOSINE_METHYLATION_RATE) * (1 - BISULFITE_CONVERSION_RATE) + CYTOSINE_METHYLATION_RATE);
                    }
                    else if(g.base2 == BaseUtils.C){
                    	POfG = (pHomVar * transversionRate) * ((1 - CYTOSINE_METHYLATION_RATE) * (1 - BISULFITE_CONVERSION_RATE) + CYTOSINE_METHYLATION_RATE);
                    }
                    else if(g.base2 == BaseUtils.T){
                    	POfG = (pTriStateGenotype * transversionRate) * ((1 - CYTOSINE_METHYLATION_RATE) * (1 - BISULFITE_CONVERSION_RATE) + CYTOSINE_METHYLATION_RATE) + (pHomVar * transversionRate) * (1 - CYTOSINE_METHYLATION_RATE) * BISULFITE_CONVERSION_RATE;
                    }
                    else{
                    	
                    }
                }
                else if(g.base1 == BaseUtils.T){
                	if(g.base2 == BaseUtils.A){
                		POfG = ((pHet - pTriStateGenotype) * transversionRate) * (1 + (1 - CYTOSINE_METHYLATION_RATE) * BISULFITE_CONVERSION_RATE);
                    }
                    else if(g.base2 == BaseUtils.G){
                    	POfG = pTriStateGenotype * transitionRate;
                    }
                    else if(g.base2 == BaseUtils.C){
                    	POfG = (pTriStateGenotype * transversionRate) * ((1 - CYTOSINE_METHYLATION_RATE) * (1 - BISULFITE_CONVERSION_RATE) + CYTOSINE_METHYLATION_RATE) + (pHomVar * transversionRate) * (1 - CYTOSINE_METHYLATION_RATE) * BISULFITE_CONVERSION_RATE;
                    }
                    else if(g.base2 == BaseUtils.T){
                    	POfG = pHomVar * transversionRate + (pTriStateGenotype * transversionRate) * (1 - CYTOSINE_METHYLATION_RATE) * BISULFITE_CONVERSION_RATE;
                    }
                    else{
                    	
                    }
                }
                else{
                	
                }
            }
            else if(ref == BaseUtils.G){
            	if(g.base1 == BaseUtils.A){
            		if(g.base2 == BaseUtils.A){
            			POfG = pHomVar * transitionRate;
                    }
                    else if(g.base2 == BaseUtils.G){
                    	POfG = (pHet - pTriStateGenotype) * transitionRate;
                    }
                    else if(g.base2 == BaseUtils.C){
                    	POfG = (pTriStateGenotype * transitionRate) * ((1 - CYTOSINE_METHYLATION_RATE) * (1 - BISULFITE_CONVERSION_RATE) + CYTOSINE_METHYLATION_RATE);
                    }
                    else if(g.base2 == BaseUtils.T){
                    	POfG = pTriStateGenotype * transitionRate;
                    }
                    else{
                    	
                    }
                }
                else if(g.base1 == BaseUtils.G){
                	if(g.base2 == BaseUtils.A){
                		POfG = (pHet - pTriStateGenotype) * transitionRate;
                    }
                    else if(g.base2 == BaseUtils.G){
                    	POfG = pHomRef;
                    }
                    else if(g.base2 == BaseUtils.C){
                    	POfG = ((pHet - pTriStateGenotype) * transversionRate) * ((1 - CYTOSINE_METHYLATION_RATE) * (1 - BISULFITE_CONVERSION_RATE) + CYTOSINE_METHYLATION_RATE);
                    }
                    else if(g.base2 == BaseUtils.T){
                    	POfG = ((pHet - pTriStateGenotype) * transversionRate) * (1 + (1 - CYTOSINE_METHYLATION_RATE) * BISULFITE_CONVERSION_RATE);
                    }
                    else{
                    	
                    }
                }
                else if(g.base1 == BaseUtils.C){
                	if(g.base2 == BaseUtils.A){
                		POfG = (pTriStateGenotype * transitionRate) * ((1 - CYTOSINE_METHYLATION_RATE) * (1 - BISULFITE_CONVERSION_RATE) + CYTOSINE_METHYLATION_RATE);
                    }
                    else if(g.base2 == BaseUtils.G){
                    	POfG = ((pHet - pTriStateGenotype) * transversionRate) * ((1 - CYTOSINE_METHYLATION_RATE) * (1 - BISULFITE_CONVERSION_RATE) + CYTOSINE_METHYLATION_RATE);
                    }
                    else if(g.base2 == BaseUtils.C){
                    	POfG = (pHomVar * transversionRate) * ((1 - CYTOSINE_METHYLATION_RATE) * (1 - BISULFITE_CONVERSION_RATE) + CYTOSINE_METHYLATION_RATE);
                    }
                    else if(g.base2 == BaseUtils.T){
                    	POfG = (pTriStateGenotype * transversionRate) * ((1 - CYTOSINE_METHYLATION_RATE) * (1 - BISULFITE_CONVERSION_RATE) + CYTOSINE_METHYLATION_RATE) + (pHomVar * transversionRate) * (1 - CYTOSINE_METHYLATION_RATE) * BISULFITE_CONVERSION_RATE;
                    }
                    else{
                    	
                    }
                }
                else if(g.base1 == BaseUtils.T){
                	if(g.base2 == BaseUtils.A){
                		POfG = pTriStateGenotype * transitionRate;
                    }
                    else if(g.base2 == BaseUtils.G){
                    	POfG = ((pHet - pTriStateGenotype) * transversionRate) * (1 + (1 - CYTOSINE_METHYLATION_RATE) * BISULFITE_CONVERSION_RATE);
                    }
                    else if(g.base2 == BaseUtils.C){
                    	POfG = (pTriStateGenotype * transversionRate) * ((1 - CYTOSINE_METHYLATION_RATE) * (1 - BISULFITE_CONVERSION_RATE) + CYTOSINE_METHYLATION_RATE) + (pHomVar * transversionRate) * (1 - CYTOSINE_METHYLATION_RATE) * BISULFITE_CONVERSION_RATE;
                    }
                    else if(g.base2 == BaseUtils.T){
                    	POfG = pHomVar * transversionRate + (pTriStateGenotype * transversionRate) * (1 - CYTOSINE_METHYLATION_RATE) * BISULFITE_CONVERSION_RATE;
                    }
                    else{
                    	
                    }
                }
                else{
                	
                }
            }
            else if(ref == BaseUtils.C){
            	if(g.base1 == BaseUtils.A){
            		if(g.base2 == BaseUtils.A){
            			POfG = pHomVar * transversionRate;
                    }
                    else if(g.base2 == BaseUtils.G){
                    	POfG = pTriStateGenotype * transversionRate;
                    }
                    else if(g.base2 == BaseUtils.C){
                    	POfG = ((pHet - pTriStateGenotype) * transversionRate) * ((1 - CYTOSINE_METHYLATION_RATE) * (1 - BISULFITE_CONVERSION_RATE) + CYTOSINE_METHYLATION_RATE);
                    }
                    else if(g.base2 == BaseUtils.T){
                    	POfG = pTriStateGenotype * transitionRate + ((pHet - pTriStateGenotype) * transitionRate) * (1 - CYTOSINE_METHYLATION_RATE) * BISULFITE_CONVERSION_RATE;
                    }
                    else{
                    	
                    }
                }
                else if(g.base1 == BaseUtils.G){
                	if(g.base2 == BaseUtils.A){
                		POfG = pTriStateGenotype * transversionRate;
                    }
                    else if(g.base2 == BaseUtils.G){
                    	POfG = pHomVar * transversionRate;
                    }
                    else if(g.base2 == BaseUtils.C){
                    	POfG = ((pHet - pTriStateGenotype) * transversionRate) * ((1 - CYTOSINE_METHYLATION_RATE) * (1 - BISULFITE_CONVERSION_RATE) + CYTOSINE_METHYLATION_RATE);
                    }
                    else if(g.base2 == BaseUtils.T){
                    	POfG = pTriStateGenotype * transitionRate + ((pHet - pTriStateGenotype) * transitionRate) * (1 - CYTOSINE_METHYLATION_RATE) * BISULFITE_CONVERSION_RATE;
                    }
                    else{
                    	
                    }
                }
                else if(g.base1 == BaseUtils.C){
                	if(g.base2 == BaseUtils.A){
                		POfG = ((pHet - pTriStateGenotype) * transversionRate) * ((1 - CYTOSINE_METHYLATION_RATE) * (1 - BISULFITE_CONVERSION_RATE) + CYTOSINE_METHYLATION_RATE);
                    }
                    else if(g.base2 == BaseUtils.G){
                    	POfG = ((pHet - pTriStateGenotype) * transversionRate) * ((1 - CYTOSINE_METHYLATION_RATE) * (1 - BISULFITE_CONVERSION_RATE) + CYTOSINE_METHYLATION_RATE);
                    }
                    else if(g.base2 == BaseUtils.C){
                    	POfG = pHomRef * ((1 - CYTOSINE_METHYLATION_RATE) * (1 - BISULFITE_CONVERSION_RATE) + CYTOSINE_METHYLATION_RATE);
                    }
                    else if(g.base2 == BaseUtils.T){
                    	POfG = ((pHet - pTriStateGenotype) * transitionRate) * ((1 - CYTOSINE_METHYLATION_RATE) * (1 - BISULFITE_CONVERSION_RATE) + CYTOSINE_METHYLATION_RATE) + pHomRef * (1 - CYTOSINE_METHYLATION_RATE) * BISULFITE_CONVERSION_RATE;
                    }
                    else{
                    	
                    }
                }
                else if(g.base1 == BaseUtils.T){
                	if(g.base2 == BaseUtils.A){
                		POfG = pTriStateGenotype * transitionRate + ((pHet - pTriStateGenotype) * transitionRate) * (1 - CYTOSINE_METHYLATION_RATE) * BISULFITE_CONVERSION_RATE;
                    }
                    else if(g.base2 == BaseUtils.G){
                    	POfG = pTriStateGenotype * transitionRate + ((pHet - pTriStateGenotype) * transitionRate) * (1 - CYTOSINE_METHYLATION_RATE) * BISULFITE_CONVERSION_RATE;
                    }
                    else if(g.base2 == BaseUtils.C){
                    	POfG = ((pHet - pTriStateGenotype) * transitionRate) * ((1 - CYTOSINE_METHYLATION_RATE) * (1 - BISULFITE_CONVERSION_RATE) + CYTOSINE_METHYLATION_RATE) + pHomRef * (1 - CYTOSINE_METHYLATION_RATE) * BISULFITE_CONVERSION_RATE;
                    }
                    else if(g.base2 == BaseUtils.T){
                    	POfG = pHomVar * transitionRate + ((pHet - pTriStateGenotype) * transitionRate) * (1 - CYTOSINE_METHYLATION_RATE) * BISULFITE_CONVERSION_RATE;
                    }
                    else{
                    	
                    }
                }
                else{
                	
                }
            }
            else if(ref == BaseUtils.T){
            	if(g.base1 == BaseUtils.A){
            		if(g.base2 == BaseUtils.A){
            			POfG = pHomVar * transversionRate;
                    }
                    else if(g.base2 == BaseUtils.G){
                    	POfG = pTriStateGenotype * transversionRate;
                    }
                    else if(g.base2 == BaseUtils.C){
                    	POfG = (pTriStateGenotype * transitionRate) * ((1 - CYTOSINE_METHYLATION_RATE) * (1 - BISULFITE_CONVERSION_RATE) + CYTOSINE_METHYLATION_RATE);
                    }
                    else if(g.base2 == BaseUtils.T){
                    	POfG = (pHet - pTriStateGenotype) * transversionRate + (pTriStateGenotype * transversionRate) * (1 - CYTOSINE_METHYLATION_RATE) * BISULFITE_CONVERSION_RATE;
                    }
                    else{
                    	
                    }
                }
                else if(g.base1 == BaseUtils.G){
                	if(g.base2 == BaseUtils.A){
                		POfG = pTriStateGenotype * transversionRate;
                    }
                    else if(g.base2 == BaseUtils.G){
                    	POfG = pHomVar * transversionRate;
                    }
                    else if(g.base2 == BaseUtils.C){
                    	POfG = (pTriStateGenotype * transitionRate) * ((1 - CYTOSINE_METHYLATION_RATE) * (1 - BISULFITE_CONVERSION_RATE) + CYTOSINE_METHYLATION_RATE);
                    }
                    else if(g.base2 == BaseUtils.T){
                    	POfG = (pHet - pTriStateGenotype) * transversionRate + (pTriStateGenotype * transversionRate) * (1 - CYTOSINE_METHYLATION_RATE) * BISULFITE_CONVERSION_RATE;
                    }
                    else{
                    	
                    }
                }
                else if(g.base1 == BaseUtils.C){
                	if(g.base2 == BaseUtils.A){
                		POfG = (pTriStateGenotype * transitionRate) * ((1 - CYTOSINE_METHYLATION_RATE) * (1 - BISULFITE_CONVERSION_RATE) + CYTOSINE_METHYLATION_RATE);
                    }
                    else if(g.base2 == BaseUtils.G){
                    	POfG = (pTriStateGenotype * transitionRate) * ((1 - CYTOSINE_METHYLATION_RATE) * (1 - BISULFITE_CONVERSION_RATE) + CYTOSINE_METHYLATION_RATE);
                    }
                    else if(g.base2 == BaseUtils.C){
                    	POfG = (pHomVar * transitionRate) * ((1 - CYTOSINE_METHYLATION_RATE) * (1 - BISULFITE_CONVERSION_RATE) + CYTOSINE_METHYLATION_RATE);
                    }
                    else if(g.base2 == BaseUtils.T){
                    	POfG = ((pHet - pTriStateGenotype) * transitionRate) * ((1 - CYTOSINE_METHYLATION_RATE) * (1 - BISULFITE_CONVERSION_RATE) + CYTOSINE_METHYLATION_RATE) + (pHomVar * transitionRate) * (1 - CYTOSINE_METHYLATION_RATE) * BISULFITE_CONVERSION_RATE;
                    }
                    else{
                    	
                    }
                }
                else if(g.base1 == BaseUtils.T){
                	if(g.base2 == BaseUtils.A){
                		POfG = (pHet - pTriStateGenotype) * transversionRate + (pTriStateGenotype * transversionRate) * (1 - CYTOSINE_METHYLATION_RATE) * BISULFITE_CONVERSION_RATE;
                    }
                    else if(g.base2 == BaseUtils.G){
                    	POfG = (pHet - pTriStateGenotype) * transversionRate + (pTriStateGenotype * transversionRate) * (1 - CYTOSINE_METHYLATION_RATE) * BISULFITE_CONVERSION_RATE;
                    }
                    else if(g.base2 == BaseUtils.C){
                    	POfG = ((pHet - pTriStateGenotype) * transitionRate) * ((1 - CYTOSINE_METHYLATION_RATE) * (1 - BISULFITE_CONVERSION_RATE) + CYTOSINE_METHYLATION_RATE) + (pHomVar * transitionRate) * (1 - CYTOSINE_METHYLATION_RATE) * BISULFITE_CONVERSION_RATE;
                    }
                    else if(g.base2 == BaseUtils.T){
                    	POfG = pHomRef + ((pHet - pTriStateGenotype) * transitionRate) * (1 - CYTOSINE_METHYLATION_RATE) * BISULFITE_CONVERSION_RATE;
                    }
                    else{
                    	
                    }
                }
                else{
                	
                }
            }
            else{
            	
            }
            
/*
            switch(ref){
    			case BaseUtils.A:{
    				switch(g.base1){
    					case BaseUtils.A:{
    						switch(g.base2){
        						//case BaseUtils.A: POfG = pHomRef;break;
        						case BaseUtils.G: POfG = (pHet - pTriStateGenotype) * transitionRate;break;
        						case BaseUtils.C: POfG = ((pHet - pTriStateGenotype) * transversionRate) * ((1 - CYTOSINE_METHYLATION_RATE) * (1 - BISULFITE_CONVERSION_RATE) + CYTOSINE_METHYLATION_RATE);break;
        						case BaseUtils.T: POfG = ((pHet - pTriStateGenotype) * transversionRate) * (1 + (1 - CYTOSINE_METHYLATION_RATE) * BISULFITE_CONVERSION_RATE);break;
        						default:POfG = pHomRef;break;
    						}
    					}break;
    					case BaseUtils.G:{
    						switch(g.base2){
        						case BaseUtils.A: POfG = (pHet - pTriStateGenotype) * transitionRate;break;
        						case BaseUtils.C: POfG = (pTriStateGenotype * transitionRate) * ((1 - CYTOSINE_METHYLATION_RATE) * (1 - BISULFITE_CONVERSION_RATE) + CYTOSINE_METHYLATION_RATE);break;
        						case BaseUtils.T: POfG = pTriStateGenotype * transitionRate;break;
        						default:POfG = pHomVar * transitionRate;break;
    						}
    					}break;
    					case BaseUtils.C:{
    						switch(g.base2){
        						case BaseUtils.A: POfG = ((pHet - pTriStateGenotype) * transversionRate) * ((1 - CYTOSINE_METHYLATION_RATE) * (1 - BISULFITE_CONVERSION_RATE) + CYTOSINE_METHYLATION_RATE);break;
        						case BaseUtils.G: POfG = (pTriStateGenotype * transitionRate) * ((1 - CYTOSINE_METHYLATION_RATE) * (1 - BISULFITE_CONVERSION_RATE) + CYTOSINE_METHYLATION_RATE);break;
        						case BaseUtils.T: POfG = (pTriStateGenotype * transversionRate) * ((1 - CYTOSINE_METHYLATION_RATE) * (1 - BISULFITE_CONVERSION_RATE) + CYTOSINE_METHYLATION_RATE) + (pHomVar * transversionRate) * (1 - CYTOSINE_METHYLATION_RATE) * BISULFITE_CONVERSION_RATE;break;
        						default:POfG = (pHomVar * transversionRate) * ((1 - CYTOSINE_METHYLATION_RATE) * (1 - BISULFITE_CONVERSION_RATE) + CYTOSINE_METHYLATION_RATE);break;
    						}
    					}break;
    					case BaseUtils.T:{
    						switch(g.base2){
        						case BaseUtils.A: POfG = ((pHet - pTriStateGenotype) * transversionRate) * (1 + (1 - CYTOSINE_METHYLATION_RATE) * BISULFITE_CONVERSION_RATE);break;
        						case BaseUtils.G: POfG = pTriStateGenotype * transitionRate;break;
        						case BaseUtils.C: POfG = (pTriStateGenotype * transversionRate) * ((1 - CYTOSINE_METHYLATION_RATE) * (1 - BISULFITE_CONVERSION_RATE) + CYTOSINE_METHYLATION_RATE) + (pHomVar * transversionRate) * (1 - CYTOSINE_METHYLATION_RATE) * BISULFITE_CONVERSION_RATE;break;
        						default:POfG = pHomVar * transversionRate + (pTriStateGenotype * transversionRate) * (1 - CYTOSINE_METHYLATION_RATE) * BISULFITE_CONVERSION_RATE;break;
    						}
    					}break;
    					default: break;
    				}
    			}break;
    			case BaseUtils.G:{
    				switch(g.base1){
    					case BaseUtils.G:{
    						switch(g.base2){
        						//case BaseUtils.A: POfG = pHomRef;break;
        						case BaseUtils.A: POfG = (pHet - pTriStateGenotype) * transitionRate;break;
        						case BaseUtils.C: POfG = ((pHet - pTriStateGenotype) * transversionRate) * ((1 - CYTOSINE_METHYLATION_RATE) * (1 - BISULFITE_CONVERSION_RATE) + CYTOSINE_METHYLATION_RATE);break;
        						case BaseUtils.T: POfG = ((pHet - pTriStateGenotype) * transversionRate) * (1 + (1 - CYTOSINE_METHYLATION_RATE) * BISULFITE_CONVERSION_RATE);break;
        						default:POfG = pHomRef;break;
    						}
    					}break;
    					case BaseUtils.A:{
    						switch(g.base2){
        						case BaseUtils.G: POfG = (pHet - pTriStateGenotype) * transitionRate;break;
        						case BaseUtils.C: POfG = (pTriStateGenotype * transitionRate) * ((1 - CYTOSINE_METHYLATION_RATE) * (1 - BISULFITE_CONVERSION_RATE) + CYTOSINE_METHYLATION_RATE);break;
        						case BaseUtils.T: POfG = pTriStateGenotype * transitionRate;break;
        						default:POfG = pHomVar * transitionRate;break;
    						}
    					}break;
    					case BaseUtils.C:{
    						switch(g.base2){
        						case BaseUtils.G: POfG = ((pHet - pTriStateGenotype) * transversionRate) * ((1 - CYTOSINE_METHYLATION_RATE) * (1 - BISULFITE_CONVERSION_RATE) + CYTOSINE_METHYLATION_RATE);break;
        						case BaseUtils.A: POfG = (pTriStateGenotype * transitionRate) * ((1 - CYTOSINE_METHYLATION_RATE) * (1 - BISULFITE_CONVERSION_RATE) + CYTOSINE_METHYLATION_RATE);break;
        						case BaseUtils.T: POfG = (pTriStateGenotype * transversionRate) * ((1 - CYTOSINE_METHYLATION_RATE) * (1 - BISULFITE_CONVERSION_RATE) + CYTOSINE_METHYLATION_RATE) + (pHomVar * transversionRate) * (1 - CYTOSINE_METHYLATION_RATE) * BISULFITE_CONVERSION_RATE;break;
        						default:POfG = (pHomVar * transversionRate) * ((1 - CYTOSINE_METHYLATION_RATE) * (1 - BISULFITE_CONVERSION_RATE) + CYTOSINE_METHYLATION_RATE);break;
    						}
    					}break;
    					case BaseUtils.T:{
    						switch(g.base2){
        						case BaseUtils.G: POfG = ((pHet - pTriStateGenotype) * transversionRate) * (1 + (1 - CYTOSINE_METHYLATION_RATE) * BISULFITE_CONVERSION_RATE);break;
        						case BaseUtils.A: POfG = pTriStateGenotype * transitionRate;break;
        						case BaseUtils.C: POfG = (pTriStateGenotype * transversionRate) * ((1 - CYTOSINE_METHYLATION_RATE) * (1 - BISULFITE_CONVERSION_RATE) + CYTOSINE_METHYLATION_RATE) + (pHomVar * transversionRate) * (1 - CYTOSINE_METHYLATION_RATE) * BISULFITE_CONVERSION_RATE;break;
        						default:POfG = pHomVar * transversionRate + (pTriStateGenotype * transversionRate) * (1 - CYTOSINE_METHYLATION_RATE) * BISULFITE_CONVERSION_RATE;break;
    						}
    					}break;
    					default: break;
    				}
    			}break;
    			case BaseUtils.C:{
    				switch(g.base1){
    					case BaseUtils.G:{
    						switch(g.base2){
        						//case BaseUtils.A: POfG = pHomRef;break;
        						case BaseUtils.A: POfG = pTriStateGenotype * transversionRate;break;
        						case BaseUtils.C: POfG = ((pHet - pTriStateGenotype) * transversionRate) * ((1 - CYTOSINE_METHYLATION_RATE) * (1 - BISULFITE_CONVERSION_RATE) + CYTOSINE_METHYLATION_RATE);break;
        						case BaseUtils.T: POfG = pTriStateGenotype * transitionRate + ((pHet - pTriStateGenotype) * transitionRate) * (1 - CYTOSINE_METHYLATION_RATE) * BISULFITE_CONVERSION_RATE;break;
        						default:POfG = pHomVar * transversionRate;break;
    						}
    					}break;
    					case BaseUtils.A:{
    						switch(g.base2){
        						case BaseUtils.G: POfG = pTriStateGenotype * transversionRate;break;
        						case BaseUtils.C: POfG = ((pHet - pTriStateGenotype) * transversionRate) * ((1 - CYTOSINE_METHYLATION_RATE) * (1 - BISULFITE_CONVERSION_RATE) + CYTOSINE_METHYLATION_RATE);break;
        						case BaseUtils.T: POfG = pTriStateGenotype * transitionRate + ((pHet - pTriStateGenotype) * transitionRate) * (1 - CYTOSINE_METHYLATION_RATE) * BISULFITE_CONVERSION_RATE;break;
        						default:POfG = pHomVar * transversionRate;break;
    						}
    					}break;
    					case BaseUtils.C:{
    						switch(g.base2){
        						case BaseUtils.G: POfG = ((pHet - pTriStateGenotype) * transversionRate) * ((1 - CYTOSINE_METHYLATION_RATE) * (1 - BISULFITE_CONVERSION_RATE) + CYTOSINE_METHYLATION_RATE);break;
        						case BaseUtils.A: POfG = ((pHet - pTriStateGenotype) * transversionRate) * ((1 - CYTOSINE_METHYLATION_RATE) * (1 - BISULFITE_CONVERSION_RATE) + CYTOSINE_METHYLATION_RATE);break;
        						case BaseUtils.T: POfG = ((pHet - pTriStateGenotype) * transitionRate) * ((1 - CYTOSINE_METHYLATION_RATE) * (1 - BISULFITE_CONVERSION_RATE) + CYTOSINE_METHYLATION_RATE) + pHomRef * (1 - CYTOSINE_METHYLATION_RATE) * BISULFITE_CONVERSION_RATE;break;
        						default:POfG = pHomRef * ((1 - CYTOSINE_METHYLATION_RATE) * (1 - BISULFITE_CONVERSION_RATE) + CYTOSINE_METHYLATION_RATE);break;
    						}
    					}break;
    					case BaseUtils.T:{
    						switch(g.base2){
        						case BaseUtils.G: POfG = pTriStateGenotype * transitionRate + ((pHet - pTriStateGenotype) * transitionRate) * (1 - CYTOSINE_METHYLATION_RATE) * BISULFITE_CONVERSION_RATE;break;
        						case BaseUtils.A: POfG = pTriStateGenotype * transitionRate + ((pHet - pTriStateGenotype) * transitionRate) * (1 - CYTOSINE_METHYLATION_RATE) * BISULFITE_CONVERSION_RATE;break;
        						case BaseUtils.C: POfG = ((pHet - pTriStateGenotype) * transitionRate) * ((1 - CYTOSINE_METHYLATION_RATE) * (1 - BISULFITE_CONVERSION_RATE) + CYTOSINE_METHYLATION_RATE) + pHomRef * (1 - CYTOSINE_METHYLATION_RATE) * BISULFITE_CONVERSION_RATE;break;
        						default:POfG = pHomVar * transitionRate + ((pHet - pTriStateGenotype) * transitionRate) * (1 - CYTOSINE_METHYLATION_RATE) * BISULFITE_CONVERSION_RATE;break;
    						}
    					}break;
    					default: break;
    				}
    			}break;
    			
    			case BaseUtils.T:{
    				switch(g.base1){
    					case BaseUtils.G:{
    						switch(g.base2){
        						//case BaseUtils.A: POfG = pHomRef;break;
        						case BaseUtils.A: POfG = pTriStateGenotype * transversionRate;break;
        						case BaseUtils.C: POfG = (pTriStateGenotype * transitionRate) * ((1 - CYTOSINE_METHYLATION_RATE) * (1 - BISULFITE_CONVERSION_RATE) + CYTOSINE_METHYLATION_RATE);break;
        						case BaseUtils.T: POfG = (pHet - pTriStateGenotype) * transversionRate + (pTriStateGenotype * transversionRate) * (1 - CYTOSINE_METHYLATION_RATE) * BISULFITE_CONVERSION_RATE;break;
        						default:POfG = pHomVar * transversionRate;break;
    						}
    					}break;
    					case BaseUtils.A:{
    						switch(g.base2){
        						case BaseUtils.G: POfG = pTriStateGenotype * transversionRate;break;
        						case BaseUtils.C: POfG = (pTriStateGenotype * transitionRate) * ((1 - CYTOSINE_METHYLATION_RATE) * (1 - BISULFITE_CONVERSION_RATE) + CYTOSINE_METHYLATION_RATE);break;
        						case BaseUtils.T: POfG = (pHet - pTriStateGenotype) * transversionRate + (pTriStateGenotype * transversionRate) * (1 - CYTOSINE_METHYLATION_RATE) * BISULFITE_CONVERSION_RATE;break;
        						default:POfG = pHomVar * transversionRate;break;
    						}
    					}break;
    					case BaseUtils.C:{
    						switch(g.base2){
        						case BaseUtils.G: POfG = (pTriStateGenotype * transitionRate) * ((1 - CYTOSINE_METHYLATION_RATE) * (1 - BISULFITE_CONVERSION_RATE) + CYTOSINE_METHYLATION_RATE);break;
        						case BaseUtils.A: POfG = (pTriStateGenotype * transitionRate) * ((1 - CYTOSINE_METHYLATION_RATE) * (1 - BISULFITE_CONVERSION_RATE) + CYTOSINE_METHYLATION_RATE);break;
        						case BaseUtils.T: POfG = ((pHet - pTriStateGenotype) * transitionRate) * ((1 - CYTOSINE_METHYLATION_RATE) * (1 - BISULFITE_CONVERSION_RATE) + CYTOSINE_METHYLATION_RATE) + (pHomVar * transitionRate) * (1 - CYTOSINE_METHYLATION_RATE) * BISULFITE_CONVERSION_RATE;break;
        						default:POfG = (pHomVar * transitionRate) * ((1 - CYTOSINE_METHYLATION_RATE) * (1 - BISULFITE_CONVERSION_RATE) + CYTOSINE_METHYLATION_RATE);break;
    						}
    					}break;
    					case BaseUtils.T:{
    						switch(g.base2){
        						case BaseUtils.G: POfG = (pHet - pTriStateGenotype) * transversionRate + (pTriStateGenotype * transversionRate) * (1 - CYTOSINE_METHYLATION_RATE) * BISULFITE_CONVERSION_RATE;break;
        						case BaseUtils.A: POfG = (pHet - pTriStateGenotype) * transversionRate + (pTriStateGenotype * transversionRate) * (1 - CYTOSINE_METHYLATION_RATE) * BISULFITE_CONVERSION_RATE;break;
        						case BaseUtils.C: POfG = ((pHet - pTriStateGenotype) * transitionRate) * ((1 - CYTOSINE_METHYLATION_RATE) * (1 - BISULFITE_CONVERSION_RATE) + CYTOSINE_METHYLATION_RATE) + (pHomVar * transitionRate) * (1 - CYTOSINE_METHYLATION_RATE) * BISULFITE_CONVERSION_RATE;break;
        						default:POfG = pHomRef + ((pHet - pTriStateGenotype) * transitionRate) * (1 - CYTOSINE_METHYLATION_RATE) * BISULFITE_CONVERSION_RATE;break;
    						}
    					}break;
    					default: break;
    				}
    			}break;
    			default: break;
            }
*/
            priors[g.ordinal()] = Math.log10(POfG);
        }

        return priors;
    }
    
   /*
    public static double[] getReferencePolarizedPriorsBasedOnMethyStatus(byte ref, double heterozyosity, double pRefError) {
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

            final double transitionRate = 2.0/3.0;
            final double transversionRate = 1.0/6.0;
            //final double nHomVars = 3;

            if ( g.isHomRef(ref) ){  
            	if(ref == BaseUtils.C){
            		//if(refNextBase == BaseUtils.G){
            			POfG = pHomRef * ((1 - CYTOSINE_METHYLATION_RATE) * (1-BISULFITE_CONVERSION_RATE) + CYTOSINE_METHYLATION_RATE);
            		//}
            		
            		
            	}
            	else if(ref == BaseUtils.T){
            		
            			POfG = pHomRef + (pHet - pTriStateGenotype) * (1 - CYTOSINE_METHYLATION_RATE) * BISULFITE_CONVERSION_RATE * transitionRate;
            		
            	}
            	else{
            		POfG = pHomRef;
            	}
            }
            else if ( g.isHomVar(ref) ){ 
            	if(ref == BaseUtils.C && g.base1 == BaseUtils.T){
            		
            			POfG = pHomVar * transitionRate + (pHet - pTriStateGenotype) * (1 - CYTOSINE_METHYLATION_RATE) * BISULFITE_CONVERSION_RATE * transitionRate;
            		
            		
            	}
            	else if((ref == BaseUtils.A || ref == BaseUtils.G || ref == BaseUtils.T) && g.base1 == BaseUtils.C){
            		
            		if(isTransition(g.base1,ref)){
            			
            				POfG = (pHomVar * transitionRate) * ((1 - CYTOSINE_METHYLATION_RATE) * (1-BISULFITE_CONVERSION_RATE) + CYTOSINE_METHYLATION_RATE);
            			
            			
            		}
            		else{
            			
            				POfG = (pHomVar * transversionRate) * ((1 - CYTOSINE_METHYLATION_RATE) * (1-BISULFITE_CONVERSION_RATE) + CYTOSINE_METHYLATION_RATE);
            			
            			
            		}
            		
            	}
            	else if((ref == BaseUtils.A || ref == BaseUtils.G) && g.base1 == BaseUtils.T){
            		
            			POfG = pHomVar * transversionRate + (pTriStateGenotype * transversionRate) * (1 - CYTOSINE_METHYLATION_RATE) * BISULFITE_CONVERSION_RATE;
            		
            		 
            	}
            	else{
            		if(isTransition(g.base1,ref)){
            			POfG = pHomVar * transitionRate;
            		}
            		else{
            			POfG = pHomVar * transversionRate;
            		}
            		
            	}
            }
            else if ( g.isHetRef(ref) ){ 
            	if(ref == BaseUtils.C){
            	     if((g.base1 == BaseUtils.C && g.base2 == BaseUtils.T)|| (g.base1 == BaseUtils.T && g.base2 == BaseUtils.C)){
            	    	 
            	    		 POfG = ((pHet - pTriStateGenotype) * transitionRate) * ((1 - CYTOSINE_METHYLATION_RATE) * (1-BISULFITE_CONVERSION_RATE) + CYTOSINE_METHYLATION_RATE) + pHomRef * (1 - CYTOSINE_METHYLATION_RATE) * BISULFITE_CONVERSION_RATE;
            	    	 
            	    	 
            	     }
            	     else{
            	    	 
            	    		 POfG = ((pHet - pTriStateGenotype) * transversionRate) * ((1 - CYTOSINE_METHYLATION_RATE) * (1-BISULFITE_CONVERSION_RATE) + CYTOSINE_METHYLATION_RATE);
            	    	 
            	    	 
            	     }
            	}
            	else if(ref == BaseUtils.T){
            		if((g.base1 == BaseUtils.C && g.base2 == BaseUtils.T)|| (g.base1 == BaseUtils.T && g.base2 == BaseUtils.C)){
            			
            				POfG = ((pHet - pTriStateGenotype) * transitionRate) * ((1 - CYTOSINE_METHYLATION_RATE) * (1-BISULFITE_CONVERSION_RATE) + CYTOSINE_METHYLATION_RATE) + (pHomVar * transitionRate) * (1 - CYTOSINE_METHYLATION_RATE) * BISULFITE_CONVERSION_RATE;
           	    	 	
            			
           	     	}
           	     	else{
           	     		
           	     			POfG = (pHet - pTriStateGenotype) * transversionRate + (pTriStateGenotype * transversionRate) * (1 - CYTOSINE_METHYLATION_RATE) * BISULFITE_CONVERSION_RATE;
           	     		
           	     		
           	     	}
            	}
            	else{
            	     if((g.base1 == ref && g.base2 == BaseUtils.C)|| (g.base1 == BaseUtils.C && g.base2 == ref)){
            	    	 
            	    		 POfG = ((pHet - pTriStateGenotype ) * transversionRate) * ((1 - CYTOSINE_METHYLATION_RATE) * (1-BISULFITE_CONVERSION_RATE) + CYTOSINE_METHYLATION_RATE);
            	    	 
            	    	 
            	     }
            	     else if((g.base1 == ref && g.base2 == BaseUtils.T)|| (g.base1 == BaseUtils.T && g.base2 == ref)){
            	    	 
            	    		 POfG = ((pHet - pTriStateGenotype ) * transversionRate) * (1 + (1 - CYTOSINE_METHYLATION_RATE) * BISULFITE_CONVERSION_RATE);
            	    	 
            	    	 
            	     }
            	     else{
            	    	 POfG = (pHet - pTriStateGenotype ) * transitionRate;
            	     }
            	}
            }
            else{ 
            	if(ref == BaseUtils.C){
            	     if((g.base1 == BaseUtils.A && g.base2 == BaseUtils.G)|| (g.base1 == BaseUtils.G && g.base2 == BaseUtils.A)){
            	    	 POfG = pTriStateGenotype * transversionRate;
            	     }
            	     else{
            	    	 
            	    		 POfG = pTriStateGenotype * transitionRate + ((pHet - pTriStateGenotype) * transitionRate) * (1 - CYTOSINE_METHYLATION_RATE) * BISULFITE_CONVERSION_RATE;
            	    	 
            	    	 
            	     }
            	}
            	else if(ref == BaseUtils.T){
            	     if((g.base1 == BaseUtils.A && g.base2 == BaseUtils.G)|| (g.base1 == BaseUtils.G && g.base2 == BaseUtils.A)){
            	    	 POfG = pTriStateGenotype * transversionRate;
            	     }
            	     else{
            	    	 
            	    		 POfG = (pTriStateGenotype * transitionRate) * ((1 - CYTOSINE_METHYLATION_RATE) * (1-BISULFITE_CONVERSION_RATE) + CYTOSINE_METHYLATION_RATE);
            	    	 
            	    	 
            	     } 
            	}
            	else{
            	     if((g.base1 == BaseUtils.C && g.base2 == BaseUtils.G)|| (g.base1 == BaseUtils.G && g.base2 == BaseUtils.C)){
            	    	 if(isTransition(g.base1,ref) || isTransition(g.base2,ref)){
            	    		 
            	    			 POfG =  (pTriStateGenotype * transitionRate) * ((1 - CYTOSINE_METHYLATION_RATE) * (1-BISULFITE_CONVERSION_RATE) + CYTOSINE_METHYLATION_RATE);
                	    	 
            	    		 
            	    	 }
            	    	 else{
            	    		 
            	    			 POfG =  (pTriStateGenotype * transversionRate) * ((1 - CYTOSINE_METHYLATION_RATE) * (1-BISULFITE_CONVERSION_RATE) + CYTOSINE_METHYLATION_RATE);
                	    	 
            	    		 
            	    	 }
            	    	 
            	     }
            	     else if((g.base1 == BaseUtils.G && g.base2 == BaseUtils.T)|| (g.base1 == BaseUtils.T && g.base2 == BaseUtils.G)){
            	    	 if(isTransition(g.base1,ref) || isTransition(g.base2,ref)){
            	    		 POfG = pTriStateGenotype * transitionRate;
            	    	 }
            	    	 else{
            	    		 POfG = pTriStateGenotype * transversionRate;
            	    	 }
            	    	 
            	     }
            	     else{
            	    	 if(isTransition(g.base1,ref) || isTransition(g.base2,ref)){
            	    		
            	    			 POfG =  (pTriStateGenotype * transitionRate) * ((1 - CYTOSINE_METHYLATION_RATE) * (1-BISULFITE_CONVERSION_RATE) + CYTOSINE_METHYLATION_RATE) + (pHomVar * transitionRate) * (1 - CYTOSINE_METHYLATION_RATE) * BISULFITE_CONVERSION_RATE;
                	    	 
            	    		 
            	    	 }
            	    	 else{
            	    		 
            	    			 POfG =  (pTriStateGenotype * transversionRate) * ((1 - CYTOSINE_METHYLATION_RATE) * (1-BISULFITE_CONVERSION_RATE) + CYTOSINE_METHYLATION_RATE) + (pHomVar * transversionRate) * (1 - CYTOSINE_METHYLATION_RATE) * BISULFITE_CONVERSION_RATE;
                	    	 
            	    		 
            	    	 }
            	    	 
            	     }
            	}
            	 
            }

            priors[g.ordinal()] = Math.log10(POfG);
        }

        return priors;
    }
*/
    static {
        for ( DiploidGenotype g : DiploidGenotype.values() ) {
            flatPriors[g.ordinal()] = Math.log10(1.0 / DiploidGenotype.values().length);
        }
    }

}
