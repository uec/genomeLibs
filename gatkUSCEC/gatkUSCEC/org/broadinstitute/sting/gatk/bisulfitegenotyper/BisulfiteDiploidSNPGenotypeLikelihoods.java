package org.broadinstitute.sting.gatk.bisulfitegenotyper;

import static java.lang.Math.log10;
import static java.lang.Math.pow;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Set;

import net.sf.samtools.SAMUtils;

import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.genotyper.DiploidSNPGenotypeLikelihoods;
import org.broadinstitute.sting.gatk.walkers.genotyper.DiploidSNPGenotypePriors;
import org.broadinstitute.sting.gatk.walkers.genotyper.PerFragmentPileupElement;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.QualityUtils;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.genotype.DiploidGenotype;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;

public class BisulfiteDiploidSNPGenotypeLikelihoods extends
		DiploidSNPGenotypeLikelihoods {
	protected BisulfiteDiploidSNPGenotypePriors priors = null;
    // TODO: don't calculate this each time through

    protected double Bisulfite_conversion_rate = 0;
    public static final double HUMAN_HETEROZYGOSITY = 1e-3;
    public static final double CEU_HETEROZYGOSITY = 1e-3;
    public static final double YRI_HETEROZYGOSITY = 1.0 / 850;
    public static final double PROB_OF_REFERENCE_ERROR = 1e-6;
    protected boolean VERBOSE = false;
    
    static BisulfiteDiploidSNPGenotypeLikelihoods[][][][][] CACHE = new BisulfiteDiploidSNPGenotypeLikelihoods[BaseUtils.BASES.length][QualityUtils.MAX_QUAL_SCORE+1][BaseUtils.BASES.length+1][QualityUtils.MAX_QUAL_SCORE+1][MAX_PLOIDY];
    
	public BisulfiteDiploidSNPGenotypeLikelihoods() {
		// TODO Auto-generated constructor stub
		super.priors = new DiploidSNPGenotypePriors();
        log10_PCR_error_3 = log10(DEFAULT_PCR_ERROR_RATE) - log10_3;
        this.priors = new BisulfiteDiploidSNPGenotypePriors();
        setToZeroBs();
	}

	public BisulfiteDiploidSNPGenotypeLikelihoods(
			DiploidSNPGenotypePriors priors, double PCR_error_rate) {
		super(priors, PCR_error_rate);
		// TODO Auto-generated constructor stub
        log10_PCR_error_3 = log10(PCR_error_rate) - log10_3;
        this.priors = new BisulfiteDiploidSNPGenotypePriors();
        setToZeroBs();
	}
	
	public BisulfiteDiploidSNPGenotypeLikelihoods(
			BisulfiteDiploidSNPGenotypePriors priors, double PCR_error_rate) {
		this.priors = priors;
        log10_1_minus_PCR_error = log10(1.0 - PCR_error_rate);
		// TODO Auto-generated constructor stub
        log10_PCR_error_3 = log10(PCR_error_rate) - log10_3;

        setToZeroBs();
	}
	
	public BisulfiteDiploidSNPGenotypeLikelihoods(
			RefMetaDataTracker tracker, ReferenceContext ref, BisulfiteDiploidSNPGenotypePriors priors, double PCR_error_rate) {
		this.priors = priors;
		this.priors.setPriors(tracker, ref, HUMAN_HETEROZYGOSITY, PROB_OF_REFERENCE_ERROR, Bisulfite_conversion_rate);
        log10_1_minus_PCR_error = log10(1.0 - PCR_error_rate);
		// TODO Auto-generated constructor stub
        log10_PCR_error_3 = log10(PCR_error_rate) - log10_3;

        setToZeroBs();
	}
	
	@Override
	/**
     * Updates likelihoods and posteriors to reflect the additional observations contained within the
     * read-based pileup up by calling add(observedBase, qualityScore) for each base / qual in the
     * pileup
     *
     * @param pileup                    read pileup
     * @param ignoreBadBases            should we ignore bad bases?
     * @param capBaseQualsAtMappingQual should we cap a base's quality by its read's mapping quality?
     * @return the number of good bases found in the pileup
     */
    public int add(ReadBackedPileup pileup, boolean ignoreBadBases, boolean capBaseQualsAtMappingQual) {
        int n = 0;

        // TODO-- move this outside the UG, e.g. to ReadBackedPileup
        // combine paired reads into a single fragment
        HashMap<String, PerFragmentPileupElement> fragments = new HashMap<String, PerFragmentPileupElement>();
        for ( PileupElement p : pileup ) {
        	if((p.getRead().getReadNegativeStrandFlag() && p.getBase() == BaseUtils.A) || (!p.getRead().getReadNegativeStrandFlag() && p.getBase() == BaseUtils.T))
        		continue;
        	Set<PileupElement> fragment = new HashSet<PileupElement>();
            String readName = p.getRead().getReadName();

            if ( fragments.containsKey(readName) )
                fragment.addAll(fragments.get(readName).getPileupElements());

            fragment.add(p);
            fragments.put(readName, new PerFragmentPileupElement(fragment));
        }

        // for each fragment, add to the likelihoods
        for ( PerFragmentPileupElement fragment : fragments.values() ) {
            n += add(fragment, ignoreBadBases, capBaseQualsAtMappingQual);
        }

        return n;
    }
	
	@Override
	public int add(PerFragmentPileupElement fragment, boolean ignoreBadBases, boolean capBaseQualsAtMappingQual) {
        // TODO-- Right now we assume that there are at most 2 reads per fragment.  This assumption is fine
        // TODO--   given the current state of next-gen sequencing, but may need to be fixed in the future.
        // TODO--   However, when that happens, we'll need to be a lot smarter about the caching we do here.
        byte observedBase1 = 0, observedBase2 = 0, qualityScore1 = 0, qualityScore2 = 0;
        for ( PileupElement p : fragment ) {
            if ( !usableBase(p, ignoreBadBases) )
                continue;

            byte qual = p.getQual();
            if ( qual > SAMUtils.MAX_PHRED_SCORE )
                throw new UserException.MalformedBAM(p.getRead(), String.format("the maximum allowed quality score is %d, but a quality of %d was observed in read %s.  Perhaps your BAM incorrectly encodes the quality scores in Sanger format; see http://en.wikipedia.org/wiki/FASTQ_format for more details", SAMUtils.MAX_PHRED_SCORE, qual, p.getRead().getReadName()));
            if ( capBaseQualsAtMappingQual )
                qual = (byte)Math.min((int)p.getQual(), p.getMappingQual());

            if ( qualityScore1 == 0 ) {
                observedBase1 = p.getBase();
                qualityScore1 = qual;
            } else {
                observedBase2 = p.getBase();
                qualityScore2 = qual;
            }
        }

        // abort early if we didn't see any good bases
        if ( observedBase1 == 0 && observedBase2 == 0 )
            return 0;

        // Just look up the cached result if it's available, or compute and store it
        BisulfiteDiploidSNPGenotypeLikelihoods gl;
        if ( ! inCache(observedBase1, qualityScore1, observedBase2, qualityScore2, FIXED_PLOIDY) ) {
            gl = calculateCachedGenotypeLikelihoods(observedBase1, qualityScore1, observedBase2, qualityScore2, FIXED_PLOIDY);
        } else {
            gl = getCachedGenotypeLikelihoods(observedBase1, qualityScore1, observedBase2, qualityScore2, FIXED_PLOIDY);
        }

        // for bad bases, there are no likelihoods
        if ( gl == null )
            return 0;

        double[] likelihoods = gl.getLikelihoods();

        for ( DiploidGenotype g : DiploidGenotype.values() ) {
            double likelihood = likelihoods[g.ordinal()];
            
            //if ( VERBOSE ) {
            //    System.out.printf("  L(%c | G=%s, Q=%d, S=%s) = %f / %f%n",
            //            observedBase, g, qualityScore, pow(10,likelihood) * 100, likelihood);
            //}

            log10Likelihoods[g.ordinal()] += likelihood;
            log10Posteriors[g.ordinal()] += likelihood;
        }

        return 1;
    }
	
	protected BisulfiteDiploidSNPGenotypeLikelihoods calculateCachedGenotypeLikelihoods(byte observedBase1, byte qualityScore1, byte observedBase2, byte qualityScore2, int ploidy) {
		BisulfiteDiploidSNPGenotypeLikelihoods gl = calculateGenotypeLikelihoods(observedBase1, qualityScore1, observedBase2, qualityScore2);
        setCache(CACHE, observedBase1, qualityScore1, observedBase2, qualityScore2, ploidy, gl);
        return gl;
    }
	
	 protected BisulfiteDiploidSNPGenotypeLikelihoods calculateGenotypeLikelihoods(byte observedBase1, byte qualityScore1, byte observedBase2, byte qualityScore2) {
	        double[] log10FourBaseLikelihoods = computeLog10Likelihoods(observedBase1, qualityScore1, observedBase2, qualityScore2);

	        try {

	        	BisulfiteDiploidSNPGenotypeLikelihoods gl = (BisulfiteDiploidSNPGenotypeLikelihoods)this.clone();
	            gl.setToZero();

	            //int ccIndex=0, ttIndex=0, ctIndex=0;
	            // we need to adjust for ploidy.  We take the raw p(obs | chrom) / ploidy, which is -log10(ploidy) in log space
	            for ( DiploidGenotype g : DiploidGenotype.values() ) {

	                // todo assumes ploidy is 2 -- should be generalized.  Obviously the below code can be turned into a loop
	                double p_base = 0.0;
	                p_base += pow(10, log10FourBaseLikelihoods[BaseUtils.simpleBaseToBaseIndex(g.base1)] - ploidyAdjustment);
	                p_base += pow(10, log10FourBaseLikelihoods[BaseUtils.simpleBaseToBaseIndex(g.base2)] - ploidyAdjustment);
	                double likelihood = log10(p_base);

	                gl.log10Likelihoods[g.ordinal()] += likelihood;
	                gl.log10Posteriors[g.ordinal()] += likelihood;
	            }
	              /*if(g.base1 == BaseUtils.C && g.base2 == BaseUtils.C){
	                	ccIndex = g.ordinal();
	            	}
	            	else if(g.base1 == BaseUtils.T && g.base2 == BaseUtils.T){
	            		ttIndex = g.ordinal();
	            	}
	            	else if((g.base1 == BaseUtils.C && g.base2 == BaseUtils.T) || (g.base1 == BaseUtils.T && g.base2 == BaseUtils.C)){
	            		ctIndex = g.ordinal();
	            	}
	            }
	            double tempcc = gl.log10Likelihoods[ccIndex];
	            double tempct = gl.log10Likelihoods[ctIndex];
	            gl.log10Likelihoods[ccIndex] += log10(1-Bisulfite_conversion_rate);
                gl.log10Posteriors[ccIndex] += log10(1-Bisulfite_conversion_rate);
                gl.log10Likelihoods[ttIndex] = pow(10, gl.log10Likelihoods[ttIndex]) + pow(10, tempct + log10(Bisulfite_conversion_rate));
                gl.log10Posteriors[ttIndex] = pow(10, gl.log10Likelihoods[ttIndex]) + pow(10, tempct + log10(Bisulfite_conversion_rate));
                gl.log10Likelihoods[ctIndex] = pow(10, tempcc + log10(Bisulfite_conversion_rate)) + pow(10, gl.log10Likelihoods[ctIndex] + log10(1-Bisulfite_conversion_rate));
                gl.log10Posteriors[ctIndex] = pow(10, tempcc + log10(Bisulfite_conversion_rate)) + pow(10, gl.log10Likelihoods[ctIndex] + log10(1-Bisulfite_conversion_rate));
	            */

	            if ( VERBOSE ) {
	                for ( DiploidGenotype g : DiploidGenotype.values() ) { System.out.printf("%s\t", g); }
	                System.out.println();
	                for ( DiploidGenotype g : DiploidGenotype.values() ) { System.out.printf("%.2f\t", gl.log10Likelihoods[g.ordinal()]); }
	                System.out.println();
	            }

	            return gl;

	         } catch ( CloneNotSupportedException e ) {
	             throw new RuntimeException(e);
	         }
	 }
	 
	 
	 protected BisulfiteDiploidSNPGenotypeLikelihoods getCachedGenotypeLikelihoods(byte observedBase1, byte qualityScore1, byte observedBase2, byte qualityScore2, int ploidy) {
		 BisulfiteDiploidSNPGenotypeLikelihoods gl = getCache(CACHE, observedBase1, qualityScore1, observedBase2, qualityScore2, ploidy);
	        if ( gl == null )
	            throw new RuntimeException(String.format("BUG: trying to fetch an unset cached genotype likelihood at base1=%c, qual1=%d, base2=%c, qual2=%d, ploidy=%d",
	                    observedBase1, qualityScore1, observedBase2, qualityScore2, ploidy));
	        return gl;
	    }

	 
	 protected void setCache( BisulfiteDiploidSNPGenotypeLikelihoods[][][][][] cache,
             byte observedBase1, byte qualityScore1, byte observedBase2, byte qualityScore2, int ploidy,
             BisulfiteDiploidSNPGenotypeLikelihoods val ) {
		 int i = BaseUtils.simpleBaseToBaseIndex(observedBase1);
		 int j = qualityScore1;
		 int k = qualityScore2 != 0 ? BaseUtils.simpleBaseToBaseIndex(observedBase2) : BaseUtils.BASES.length;
		 int l = qualityScore2;
		 int m = ploidy;

		 cache[i][j][k][l][m] = val;
	 }
	 
	 protected BisulfiteDiploidSNPGenotypeLikelihoods getCache(BisulfiteDiploidSNPGenotypeLikelihoods[][][][][] cache,
             byte observedBase1, byte qualityScore1, byte observedBase2, byte qualityScore2, int ploidy) {
		 int i = BaseUtils.simpleBaseToBaseIndex(observedBase1);
		 int j = qualityScore1;
		 int k = qualityScore2 != 0 ? BaseUtils.simpleBaseToBaseIndex(observedBase2) : BaseUtils.BASES.length;
		 int l = qualityScore2;
		 int m = ploidy;
		 return cache[i][j][k][l][m];
	 }

	

    @Override
    public double[] getPriors() {
        return this.priors.getPriors();
    }
    
    public void setPriors(BisulfiteDiploidSNPGenotypePriors priors) {
        this.priors = priors;
        log10Posteriors = genotypeZeros.clone();
        for ( DiploidGenotype g : DiploidGenotype.values() ) {
            int i = g.ordinal();
            log10Posteriors[i] = priors.getPriors()[i] + log10Likelihoods[i];
        }
    }
 // likelihoods are all zeros

    protected void setToZeroBs() {
    	//System.err.println("test!!!");
    	log10Likelihoods = genotypeZeros.clone();
       // if(this.priors.getPriors() == null){
        //	System.err.println("priors null!!!");
       // }
       // System.err.println("log10Likelihoods null!!!");
        log10Posteriors = this.priors.getPriors().clone();     // posteriors are all the priors
    }

    protected Object clone() throws CloneNotSupportedException {
    	BisulfiteDiploidSNPGenotypeLikelihoods c = (BisulfiteDiploidSNPGenotypeLikelihoods)super.clone();
        c.priors = priors;
        c.log10Likelihoods = log10Likelihoods.clone();
        c.log10Posteriors = log10Posteriors.clone();
        return c;
    }
    
    /**
     * Sets the priors
     * @param priors priors
     */
	
    private final static double[] genotypeZeros = new double[DiploidGenotype.values().length];
	private final static double[] baseZeros = new double[BaseUtils.BASES.length];
}
