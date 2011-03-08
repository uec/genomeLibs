package org.broadinstitute.sting.gatk.bisulfitegenotyper;

import static java.lang.Math.log10;
import static java.lang.Math.pow;

import org.broadinstitute.sting.gatk.walkers.genotyper.DiploidSNPGenotypeLikelihoods;
import org.broadinstitute.sting.gatk.walkers.genotyper.DiploidSNPGenotypePriors;
import org.broadinstitute.sting.utils.BaseUtils;

public class BisulfiteDiploidSNPGenotypeLikelihoods extends
		DiploidSNPGenotypeLikelihoods {

    // TODO: don't calculate this each time through

    protected double Bisulfite_conversion_rate = 0.5;
    
	public BisulfiteDiploidSNPGenotypeLikelihoods() {
		// TODO Auto-generated constructor stub
		super.priors = new DiploidSNPGenotypePriors();
        log10_PCR_error_3 = log10(DEFAULT_PCR_ERROR_RATE) - log10_3;

        setToZero();
	}

	public BisulfiteDiploidSNPGenotypeLikelihoods(
			DiploidSNPGenotypePriors priors, double PCR_error_rate) {
		super(priors, PCR_error_rate);
		// TODO Auto-generated constructor stub
        log10_PCR_error_3 = log10(PCR_error_rate) - log10_3;

        setToZero();
	}

	@Override
	/**
     * Updates likelihoods and posteriors to reflect an additional observation of observedBase with
     * qualityScore.
     *
     * @param observedBase1  the base observed on the 1st read of the fragment
     * @param qualityScore1  the qual of the base on the 1st read of the fragment, or zero if NA
     * @param observedBase2  the base observed on the 2nd read of the fragment
     * @param qualityScore2  the qual of the base on the 2nd read of the fragment, or zero if NA
     * @return likelihoods for this observation or null if the base was not considered good enough to add to the likelihoods (Q0 or 'N', for example)
     */
    protected double[] computeLog10Likelihoods(byte observedBase1, byte qualityScore1, byte observedBase2, byte qualityScore2) {
		double[] log10FourBaseLikelihoods = baseZeros.clone();

        for ( byte trueBase : BaseUtils.BASES ) {
            double likelihood = 0.0;

            for ( byte fragmentBase : BaseUtils.BASES ) {
            	double log10FragmentLikelihood = 0;
            	if(trueBase == fragmentBase){
            		if(trueBase == BaseUtils.C || trueBase == BaseUtils.T){
            			log10FragmentLikelihood = log10_1_minus_PCR_error + log10(1.0-Bisulfite_conversion_rate);
            		}
            		else{
            			log10FragmentLikelihood = log10_1_minus_PCR_error;
            		}
            	}
            	else{
            		if((trueBase == BaseUtils.C && fragmentBase == BaseUtils.T)||(trueBase == BaseUtils.T && fragmentBase == BaseUtils.C)){
            			log10FragmentLikelihood = log10(pow(10, log10_1_minus_PCR_error + log10(Bisulfite_conversion_rate)) + pow(10, log10_PCR_error_3));
            		}
            		else{
            			log10FragmentLikelihood = log10_PCR_error_3;
            		}
            	}

                if ( qualityScore1 != 0 ) {
                    log10FragmentLikelihood += log10PofObservingBaseGivenChromosome(observedBase1, fragmentBase, qualityScore1);
                }
                if ( qualityScore2 != 0 ) {
                    log10FragmentLikelihood += log10PofObservingBaseGivenChromosome(observedBase2, fragmentBase, qualityScore2);
                }

                //if ( VERBOSE ) {
                //    System.out.printf("  L(%c | b=%s, Q=%d) = %f / %f%n",
                //            observedBase, trueBase, qualityScore, pow(10,likelihood) * 100, likelihood);
                //}

                likelihood += pow(10, log10FragmentLikelihood);
            }

            log10FourBaseLikelihoods[BaseUtils.simpleBaseToBaseIndex(trueBase)] = log10(likelihood);
        }

        return log10FourBaseLikelihoods;
	}
	
	@Override
	/**
    *
    * @param observedBase observed base
    * @param chromBase    target base
    * @param qual         base quality
    * @return log10 likelihood
    */
   protected double log10PofObservingBaseGivenChromosome(byte observedBase, byte chromBase, byte qual) {

       double logP;
       double e = pow(10, (qual / -10.0));
       if ( observedBase == chromBase ) {
           // the base is consistent with the chromosome -- it's 1 - e
           //logP = oneMinusData[qual];
    	   if(observedBase == BaseUtils.C || observedBase == BaseUtils.T){
    		   logP = log10((1.0 - e)*(1-Bisulfite_conversion_rate));
   			}
   			else{
   				logP = log10(1.0 - e);
   			}
    	   
           
       } else {
           // the base is inconsistent with the chromosome -- it's e * P(chromBase | observedBase is an error)
    	   if((observedBase == BaseUtils.C && chromBase == BaseUtils.T)||(observedBase == BaseUtils.T && chromBase == BaseUtils.C)){
    		   logP = log10((1.0 - e)*Bisulfite_conversion_rate+e/3.0);
           }
           else{
        	   logP = qual / -10.0 + (-log10_3);
           }
    	   
       }

       //System.out.printf("%c %c %d => %f%n", observedBase, chromBase, qual, logP);
       return logP;
   }
	
	private final static double[] baseZeros = new double[BaseUtils.BASES.length];
}
