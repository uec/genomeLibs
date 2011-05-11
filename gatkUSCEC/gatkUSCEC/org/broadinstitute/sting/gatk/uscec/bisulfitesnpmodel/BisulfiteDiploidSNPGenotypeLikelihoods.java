package org.broadinstitute.sting.gatk.uscec.bisulfitesnpmodel;

import static java.lang.Math.log10;
import static java.lang.Math.pow;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Set;

import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMUtils;

import org.broadinstitute.sting.gatk.uscec.bisulfitesnpmodel.BisulfiteDiploidSNPGenotypePriors;
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

    protected double BISULFITE_CONVERSION_RATE;
    protected static double CPG_METHYLATION_RATE = 0;
    protected static double CPH_METHYLATION_RATE = 0;
    protected static String DETERMINED_CYTOSINE_TYPE = "C";
    protected static double DETERMINED_CYTOSINE_TYPE_C_METHY = 0;
    protected static boolean DETERMINED_CYTOSINE_TYPE_NEGSTRAND = false;
    
    public static final double HUMAN_HETEROZYGOSITY = 1e-3;
    public static final double CEU_HETEROZYGOSITY = 1e-3;
    public static final double YRI_HETEROZYGOSITY = 1.0 / 850;
    public static final double PROB_OF_REFERENCE_ERROR = 1e-6;
    public boolean VERBOSE = false;
    
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

        this.priors = new BisulfiteDiploidSNPGenotypePriors();
        setToZeroBs();
	}
	
	public BisulfiteDiploidSNPGenotypeLikelihoods(
			BisulfiteDiploidSNPGenotypePriors priors, double PCR_error_rate) {
		this.priors = priors;

        setToZeroBs();
	}
	
	public BisulfiteDiploidSNPGenotypeLikelihoods(
			RefMetaDataTracker tracker, ReferenceContext ref, BisulfiteDiploidSNPGenotypePriors priors, double PCR_error_rate, double bisulfiteConversionRate, double cpgMethyRate, double cphMethyRate, double novelDbsnpHet, double validateDbsnpHet, byte[] contextSeq) {
		this.priors = priors;
		this.BISULFITE_CONVERSION_RATE = bisulfiteConversionRate;
		this.CPG_METHYLATION_RATE = cpgMethyRate;
		this.CPH_METHYLATION_RATE = cphMethyRate;
		DETERMINED_CYTOSINE_TYPE_C_METHY = CPH_METHYLATION_RATE;
		this.priors.setPriors(tracker, ref, HUMAN_HETEROZYGOSITY, PROB_OF_REFERENCE_ERROR, BISULFITE_CONVERSION_RATE, CPG_METHYLATION_RATE, CPH_METHYLATION_RATE, novelDbsnpHet, validateDbsnpHet, contextSeq );
        setToZeroBs();
	}
	
	public String determinCytosineStatus(ReadBackedPileup pileup, String cytosineList){
		HashMap<String, Double> cytosineListMap = new HashMap<String, Double>();
		
		String[] tmpArray = cytosineList.split(";");
		for(String tmp : tmpArray){
			String[] mapElement = tmp.split(":");
			cytosineListMap.put(mapElement[0], Double.parseDouble(mapElement[1]));
		}
		
		for(String cytosineType : cytosineListMap.keySet()){
			String[] tmpKey = cytosineType.split("-");
			Integer cytosinePos = Integer.parseInt(tmpKey[1]);
			double[] adjacentCytosineSeqLikelihood = new double[tmpKey[0].length()];
			//double[] adjacentCytosineSeqLikelihoodReverseStrand = new double[tmpKey[0].length()];
			ReadBackedPileup pileupPositiveStrand = pileup.getPositiveStrandPileup();
	        for(PileupElement p : pileupPositiveStrand ){
	        	SAMRecord samRecord = p.getRead();
	        	
				int cytosineOffset = p.getOffset();
				for(int i = 0; i < tmpKey[0].length(); i++){
					if(i == cytosinePos)
						continue;
					int elementOffset = i - cytosinePos + cytosineOffset;
					if(elementOffset < 0 || elementOffset > samRecord.getReadLength()-1)
						continue;
					//System.err.println("1: " + elementOffset + "\t" + cytosinePos + "\t" + cytosineOffset + "\t" + samRecord.getAlignmentStart() );
					PileupElement tmpP = new PileupElement(samRecord,elementOffset);
					if ( !usableBase(tmpP, true) ){
			        	continue;
			        }
					byte qual = tmpP.getQual();
					if ( qual > SAMUtils.MAX_PHRED_SCORE )
			            throw new UserException.MalformedBam(tmpP.getRead(), String.format("the maximum allowed quality score is %d, but a quality of %d was observed in read %s.  Perhaps your BAM incorrectly encodes the quality scores in Sanger format; see http://en.wikipedia.org/wiki/FASTQ_format for more details", SAMUtils.MAX_PHRED_SCORE, qual, tmpP.getRead().getReadName()));
			        qual = (byte)Math.min((int)tmpP.getQual(), tmpP.getMappingQual());
			        byte observedBase = tmpP.getBase();
			        byte qualityScore = qual;
			        byte queryByte = (byte)tmpKey[0].charAt(i);
			        double error = pow(10, -qualityScore/10.0);
			        if ( observedBase == 0 )
			            continue;
			        adjacentCytosineSeqLikelihood[i] += observedBase == queryByte ? log10(1.0 - error) : ( (queryByte == BaseUtils.C && observedBase == BaseUtils.T) ? log10(1.0 - error) : log10(error/3.0));
			        
				}
				/*
				//int cytosineOffsetOpposite = samRecord.getReadLength() - p.getOffset();
				for(int i = 0; i < tmpKey[0].length(); i++){
					if(i == cytosinePos)
						continue;
					int elementOffset = i - cytosinePos + cytosineOffset;
					if(elementOffset < 0 || elementOffset > samRecord.getReadLength()-1)
						continue;
					//System.err.println("2: " + elementOffset + "\t" + cytosinePos + "\t" + cytosineOffset + "\t" + samRecord.getAlignmentStart() );
					PileupElement tmpP = new PileupElement(samRecord,elementOffset);
					if ( !usableBase(tmpP, true) ){
			        	continue;
			        }
					byte qual = tmpP.getQual();
					if ( qual > SAMUtils.MAX_PHRED_SCORE )
			            throw new UserException.MalformedBAM(tmpP.getRead(), String.format("the maximum allowed quality score is %d, but a quality of %d was observed in read %s.  Perhaps your BAM incorrectly encodes the quality scores in Sanger format; see http://en.wikipedia.org/wiki/FASTQ_format for more details", SAMUtils.MAX_PHRED_SCORE, qual, tmpP.getRead().getReadName()));
			        qual = (byte)Math.min((int)tmpP.getQual(), tmpP.getMappingQual());
			        byte observedBase = tmpP.getBase();
			        byte qualityScore = qual;
			        byte queryByte = (byte)tmpKey[0].charAt(tmpKey[0].length() - 1 - i);
			        queryByte = BaseUtils.simpleComplement(queryByte);
			        double error = pow(10, -qualityScore/10.0);
			        if ( observedBase == 0 )
			            continue;
			        adjacentCytosineSeqLikelihoodReverseStrand[i] += observedBase == queryByte ? log10(1.0 - error) : ( (queryByte == BaseUtils.G && observedBase == BaseUtils.A) ? log10(1.0 - error) : log10(error/3.0));
			        
				}
				
		          */  	
		        
		        
	        }
	        ReadBackedPileup pileupNegativeStrand = pileup.getNegativeStrandPileup();
	        for(PileupElement p : pileupNegativeStrand ){
	        	SAMRecord samRecord = p.getRead();
	        	
				int cytosineOffset = p.getOffset();
				for(int i = 0; i < tmpKey[0].length(); i++){
					if(i == cytosinePos)
						continue;
					int elementOffset = i - cytosinePos + cytosineOffset;
					if(elementOffset < 0 || elementOffset > samRecord.getReadLength()-1)
						continue;
					//System.err.println("3: " + elementOffset + "\t" + cytosinePos + "\t" + cytosineOffset + "\t" + samRecord.getAlignmentStart() );
					PileupElement tmpP = new PileupElement(samRecord,elementOffset);
					if ( !usableBase(tmpP, true) ){
			        	continue;
			        }
					byte qual = tmpP.getQual();
					if ( qual > SAMUtils.MAX_PHRED_SCORE )
			            throw new UserException.MalformedBam(tmpP.getRead(), String.format("the maximum allowed quality score is %d, but a quality of %d was observed in read %s.  Perhaps your BAM incorrectly encodes the quality scores in Sanger format; see http://en.wikipedia.org/wiki/FASTQ_format for more details", SAMUtils.MAX_PHRED_SCORE, qual, tmpP.getRead().getReadName()));
			        qual = (byte)Math.min((int)tmpP.getQual(), tmpP.getMappingQual());
			        byte observedBase = tmpP.getBase();
			        byte qualityScore = qual;
			        byte queryByte = (byte)tmpKey[0].charAt(i);
			        double error = pow(10, -qualityScore/10.0);
			        if ( observedBase == 0 )
			            continue;
			        adjacentCytosineSeqLikelihood[i] += observedBase == queryByte ? log10(1.0 - error) : ( (queryByte == BaseUtils.G && observedBase == BaseUtils.A) ? log10(1.0 - error) : log10(error/3.0));
			        
				}
				/*
				//int cytosineOffsetOpposite = samRecord.getReadLength() - p.getOffset();
				for(int i = 0; i < tmpKey[0].length(); i++){
					if(i == cytosinePos)
						continue;
					int elementOffset = i - cytosinePos + cytosineOffset;
					if(elementOffset < 0 || elementOffset > samRecord.getReadLength()-1)
						continue;
					//System.err.println("4: " + elementOffset + "\t" + cytosinePos + "\t" + cytosineOffset + "\t" + samRecord.getAlignmentStart() );
					PileupElement tmpP = new PileupElement(samRecord,elementOffset);
					//System.err.println("4: " + tmpP.toString() );
					if ( !usableBase(tmpP, true) ){
			        	continue;
			        }
					byte qual = tmpP.getQual();
					if ( qual > SAMUtils.MAX_PHRED_SCORE )
			            throw new UserException.MalformedBAM(tmpP.getRead(), String.format("the maximum allowed quality score is %d, but a quality of %d was observed in read %s.  Perhaps your BAM incorrectly encodes the quality scores in Sanger format; see http://en.wikipedia.org/wiki/FASTQ_format for more details", SAMUtils.MAX_PHRED_SCORE, qual, tmpP.getRead().getReadName()));
			        qual = (byte)Math.min((int)tmpP.getQual(), tmpP.getMappingQual());
			        byte observedBase = tmpP.getBase();
			        byte qualityScore = qual;
			        byte queryByte = (byte)tmpKey[0].charAt(tmpKey[0].length() - 1 - i);
			        queryByte = BaseUtils.simpleComplement(queryByte);
			        double error = pow(10, -qualityScore/10.0);
			        if ( observedBase == 0 )
			            continue;
			        adjacentCytosineSeqLikelihoodReverseStrand[i] += observedBase == queryByte ? log10(1.0 - error) : ( (queryByte == BaseUtils.C && observedBase == BaseUtils.T) ? log10(1.0 - error) : log10(error/3.0));
			        
				}
				*/
	        }
	     
	        double sum = 0;
	      //  double sumReverse = 0;
	        for(int i = 0; i < adjacentCytosineSeqLikelihood.length; i++)
	        	sum += adjacentCytosineSeqLikelihood[i];
	      //  for(int i = 0; i < adjacentCytosineSeqLikelihoodReverseStrand.length; i++)
	        //	sumReverse += adjacentCytosineSeqLikelihoodReverseStrand[i];
	      //  if(sum >= log10(0.90) || sumReverse >= log10(0.90)){
	        if(sum >= log10(0.90)){
	        	DETERMINED_CYTOSINE_TYPE = tmpKey[0];
	        	DETERMINED_CYTOSINE_TYPE_C_METHY = cytosineListMap.get(cytosineType);
	        	if(sum >=log10(0.90)){
	        		DETERMINED_CYTOSINE_TYPE_NEGSTRAND = false;
	        	}
	        	else{
	        		DETERMINED_CYTOSINE_TYPE_NEGSTRAND = true;
	        	}
	        	break;
	        }
	        
		}
		
		if(VERBOSE){
			 System.out.println("DETERMINED_CYTOSINE_TYPE: " + DETERMINED_CYTOSINE_TYPE + "\tDETERMINED_CYTOSINE_TYPE_NEGSTRAND: " + DETERMINED_CYTOSINE_TYPE_NEGSTRAND + "\tDETERMINED_CYTOSINE_TYPE_C_METHY: " + DETERMINED_CYTOSINE_TYPE_C_METHY);
		}
		
		
		return DETERMINED_CYTOSINE_TYPE;
	}
	
	
	public int add(ReadBackedPileup pileup, boolean ignoreBadBases, boolean capBaseQualsAtMappingQual, byte refNextBase, byte refPreBase) {
        int n = 0;
        
        ReadBackedPileup pileupPositiveStrand = pileup.getPositiveStrandPileup();
        for(PileupElement p : pileupPositiveStrand )
        	n += add(p, ignoreBadBases, capBaseQualsAtMappingQual, false, refNextBase, refPreBase);
        
        ReadBackedPileup pileupNegativeStrand = pileup.getNegativeStrandPileup();
        for(PileupElement p : pileupNegativeStrand )
        	n += add(p, ignoreBadBases, capBaseQualsAtMappingQual, true, refNextBase, refPreBase);
 
        return n;
    }
	
	@Override
	public int add(ReadBackedPileup pileup, boolean ignoreBadBases, boolean capBaseQualsAtMappingQual) {
        int n = 0;
        
        ReadBackedPileup pileupPositiveStrand = pileup.getPositiveStrandPileup();
        for(PileupElement p : pileupPositiveStrand )
        	n += add(p, ignoreBadBases, capBaseQualsAtMappingQual, false);
        
        ReadBackedPileup pileupNegativeStrand = pileup.getNegativeStrandPileup();
        for(PileupElement p : pileupNegativeStrand )
        	n += add(p, ignoreBadBases, capBaseQualsAtMappingQual, true);
 
        return n;
    }
	
	public int add(PileupElement p, boolean ignoreBadBases, boolean capBaseQualsAtMappingQual, boolean negStrand) {
        // TODO-- Right now we assume that there are at most 2 reads per fragment.  This assumption is fine
        // TODO--   given the current state of next-gen sequencing, but may need to be fixed in the future.
        // TODO--   However, when that happens, we'll need to be a lot smarter about the caching we do here.
        byte observedBase = 0, qualityScore = 0;
        if ( !usableBase(p, ignoreBadBases) ){
        	//System.out.println("bad base: " + p.getBase());
        	return 0;
        }
            	
        byte qual = p.getQual();
        if ( qual > SAMUtils.MAX_PHRED_SCORE )
            throw new UserException.MalformedBam(p.getRead(), String.format("the maximum allowed quality score is %d, but a quality of %d was observed in read %s.  Perhaps your BAM incorrectly encodes the quality scores in Sanger format; see http://en.wikipedia.org/wiki/FASTQ_format for more details", SAMUtils.MAX_PHRED_SCORE, qual, p.getRead().getReadName()));
        if ( capBaseQualsAtMappingQual )
            qual = (byte)Math.min((int)p.getQual(), p.getMappingQual());
         
        observedBase = p.getBase();
        qualityScore = qual;
         
        // abort early if we didn't see any good bases
        if ( observedBase == 0 )
            return 0;


        BisulfiteDiploidSNPGenotypeLikelihoods gl;
        gl = getCalculateGenotypeLikelihoods(observedBase, qualityScore, negStrand);


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

	public int add(PileupElement p, boolean ignoreBadBases, boolean capBaseQualsAtMappingQual, boolean negStrand, byte refNextBase, byte refPreBase) {
        // TODO-- Right now we assume that there are at most 2 reads per fragment.  This assumption is fine
        // TODO--   given the current state of next-gen sequencing, but may need to be fixed in the future.
        // TODO--   However, when that happens, we'll need to be a lot smarter about the caching we do here.
        byte observedBase = 0, qualityScore = 0;
        if ( !usableBase(p, ignoreBadBases) ){
        	//System.out.println("bad base: " + p.getBase());
        	return 0;
        }
            	
        byte qual = p.getQual();
        if ( qual > SAMUtils.MAX_PHRED_SCORE )
            throw new UserException.MalformedBam(p.getRead(), String.format("the maximum allowed quality score is %d, but a quality of %d was observed in read %s.  Perhaps your BAM incorrectly encodes the quality scores in Sanger format; see http://en.wikipedia.org/wiki/FASTQ_format for more details", SAMUtils.MAX_PHRED_SCORE, qual, p.getRead().getReadName()));
        if ( capBaseQualsAtMappingQual )
            qual = (byte)Math.min((int)p.getQual(), p.getMappingQual());
         
        observedBase = p.getBase();
        qualityScore = qual;
         
        // abort early if we didn't see any good bases
        if ( observedBase == 0 )
            return 0;


        BisulfiteDiploidSNPGenotypeLikelihoods gl;
        gl = getCalculateGenotypeLikelihoods(observedBase, qualityScore, negStrand, refNextBase, refPreBase);


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
	
	/*
	 * for reads: AAGGCCTTaaggcctt (Q=30)
	 * likelihood of all of the genotype would be (err = pow(10,Q/-10))
	 * AA: (1-err)^2 * (err/3)^10 * (1-err+err*conv/3)^2 * (err/3 * (1-conv))^2
	 * AG: [1/2 * (1-err) + 1/2 * err/3]^4 *  [1/2 * err/3 + 1/2 * (1-err)]^2 * (err/3)^8 * [1/2 * (1-err) + 1/2 * (err/3) * conv]
	 * AC: [1/2 * (1-err) + 1/2 * err/3 * (1-conv)]^2 * (err/3)^2 * [1/2 * (err/3) * (1-conv) + 1/2 * (1-err)* (1-conv)]^2 * [1/2 * err/3 + 1/2 * (err/3 + conv)]^2 * [1/2 * (1-err) + 1/2 * (err/3 + conv)]^2 * (err/3)^2 * [1/2 * err/3 + 1/2 * (1-err)]^2 * (err/3)^2 
	 * AT: 
	 * 
	 * 
	*/
	protected BisulfiteDiploidSNPGenotypeLikelihoods getCalculateGenotypeLikelihoods(byte observedBase, byte qualityScore, boolean negStrand, byte refNextBase, byte refPreBase){

		try {
			BisulfiteDiploidSNPGenotypeLikelihoods gl = (BisulfiteDiploidSNPGenotypeLikelihoods)this.clone();
			gl.setToZeroBs();
			double error = pow(10, -qualityScore/10.0);
			double pOfBase1 = 0.5, pOfBase2 = 1.0 - pOfBase1;
			for ( DiploidGenotype g : DiploidGenotype.values() ) {
				double likelihood = 0.0;
				if(!negStrand){
						switch(observedBase){
							case BaseUtils.C:
								if(refNextBase == BaseUtils.G){
									likelihood += observedBase == g.base1 ? pOfBase1 * (1.0-error) * ( (1.0-CPG_METHYLATION_RATE) * (1.0-BISULFITE_CONVERSION_RATE) + CPG_METHYLATION_RATE ) : pOfBase1 * (error/3.0) * ( (1.0-CPG_METHYLATION_RATE) * (1.0-BISULFITE_CONVERSION_RATE) + CPG_METHYLATION_RATE );
									likelihood += observedBase == g.base2 ? pOfBase2 * (1.0-error) *  ( (1.0-CPG_METHYLATION_RATE) * (1.0-BISULFITE_CONVERSION_RATE) + CPG_METHYLATION_RATE ) : pOfBase2 * (error/3.0) * ( (1.0-CPG_METHYLATION_RATE) * (1.0-BISULFITE_CONVERSION_RATE) + CPG_METHYLATION_RATE );
									if ( VERBOSE ) {
									//	System.out.println("flag1: observedBase-" + observedBase + "\t" + "g.base1-" + g.base1 + "\t" + "g.base2-" + g.base2 + "\t" + "likelihood-" + log10(likelihood));
									}
								}
								else{
									likelihood += observedBase == g.base1 ? pOfBase1 * (1.0-error) * ((1.0-CPH_METHYLATION_RATE) * (1.0-BISULFITE_CONVERSION_RATE) + CPH_METHYLATION_RATE ) : pOfBase1 * (error/3.0) * ((1.0-CPH_METHYLATION_RATE) * (1.0-BISULFITE_CONVERSION_RATE) + CPH_METHYLATION_RATE );
									likelihood += observedBase == g.base2 ? pOfBase2 * (1.0-error) * ((1.0-CPH_METHYLATION_RATE) * (1.0-BISULFITE_CONVERSION_RATE) + CPH_METHYLATION_RATE ) : pOfBase2 * (error/3.0) * ((1.0-CPH_METHYLATION_RATE) * (1.0-BISULFITE_CONVERSION_RATE) + CPH_METHYLATION_RATE );
									if ( VERBOSE ) {
									//	System.out.println("flag2: observedBase-" + observedBase + "\t" + "g.base1-" + g.base1 + "\t" + "g.base2-" + g.base2 + "\t" + "likelihood-" + log10(likelihood));
									}
								}
								break;
							case BaseUtils.T:
								if(refNextBase == BaseUtils.G){
									likelihood += observedBase == g.base1 ? pOfBase1 * (1.0-error) : (g.base1 == BaseUtils.C ? pOfBase1 * (error/3.0 + (1.0-CPG_METHYLATION_RATE) * BISULFITE_CONVERSION_RATE) : pOfBase1 * (error/3.0));
									likelihood += observedBase == g.base2 ? pOfBase1 * (1.0-error) : (g.base2 == BaseUtils.C ? pOfBase2 * (error/3.0 + (1.0-CPG_METHYLATION_RATE) * BISULFITE_CONVERSION_RATE) : pOfBase2 * (error/3.0));
									if ( VERBOSE ) {
									//	System.out.println("flag3: observedBase-" + observedBase + "\t" + "g.base1-" + g.base1 + "\t" + "g.base2-" + g.base2 + "\t" + "likelihood-" + log10(likelihood));
									}
								}
								else{
									likelihood += observedBase == g.base1 ? pOfBase1 * (1.0-error) : (g.base1 == BaseUtils.C ? pOfBase1 * (error/3.0 + (1.0-CPH_METHYLATION_RATE) * BISULFITE_CONVERSION_RATE) : pOfBase1 * (error/3.0));
									likelihood += observedBase == g.base2 ? pOfBase1 * (1.0-error) : (g.base2 == BaseUtils.C ? pOfBase2 * (error/3.0 + (1.0-CPH_METHYLATION_RATE) * BISULFITE_CONVERSION_RATE) : pOfBase2 * (error/3.0));
									if ( VERBOSE ) {
									//	System.out.println("flag4: observedBase-" + observedBase + "\t" + "g.base1-" + g.base1 + "\t" + "g.base2-" + g.base2 + "\t" + "likelihood-" + log10(likelihood));
									}
								}
								
								break;
							default: 
								likelihood += observedBase == g.base1 ? pOfBase1 * (1.0-error) : pOfBase1 * (error/3.0);
								likelihood += observedBase == g.base2 ? pOfBase2 * (1.0-error) : pOfBase2 * (error/3.0);
								if ( VERBOSE ) {
									//System.out.println("flag5: observedBase-" + observedBase + "\t" + "g.base1-" + g.base1 + "\t" + "g.base2-" + g.base2 + "\t" + "likelihood-" + log10(likelihood));
								}
									
						}
				}
				else{
						switch(observedBase){
						case BaseUtils.G:
							if(refPreBase == BaseUtils.C){
								likelihood += observedBase == g.base1 ? pOfBase1 * (1.0-error) * ( (1.0-CPG_METHYLATION_RATE) * (1.0-BISULFITE_CONVERSION_RATE) + CPG_METHYLATION_RATE ) : pOfBase1 * (error/3.0) * ( (1.0-CPG_METHYLATION_RATE) * (1.0-BISULFITE_CONVERSION_RATE) + CPG_METHYLATION_RATE );
								likelihood += observedBase == g.base2 ? pOfBase2 * (1.0-error) *  ( (1.0-CPG_METHYLATION_RATE) * (1.0-BISULFITE_CONVERSION_RATE) + CPG_METHYLATION_RATE ) : pOfBase2 * (error/3.0) * ( (1.0-CPG_METHYLATION_RATE) * (1.0-BISULFITE_CONVERSION_RATE) + CPG_METHYLATION_RATE );
								if ( VERBOSE ) {
									//System.out.println("flag6: observedBase-" + observedBase + "\t" + "g.base1-" + g.base1 + "\t" + "g.base2-" + g.base2 + "\t" + "likelihood-" + log10(likelihood));
								}
							}
							else{
								likelihood += observedBase == g.base1 ? pOfBase1 * (1.0-error) * ((1.0-CPH_METHYLATION_RATE) * (1.0-BISULFITE_CONVERSION_RATE) + CPH_METHYLATION_RATE) : pOfBase1 * (error/3.0) * ((1.0-CPH_METHYLATION_RATE) * (1.0-BISULFITE_CONVERSION_RATE) + CPH_METHYLATION_RATE);
								likelihood += observedBase == g.base2 ? pOfBase2 * (1.0-error) * ((1.0-CPH_METHYLATION_RATE) * (1.0-BISULFITE_CONVERSION_RATE) + CPH_METHYLATION_RATE) : pOfBase2 * (error/3.0) * ((1.0-CPH_METHYLATION_RATE) * (1.0-BISULFITE_CONVERSION_RATE) + CPH_METHYLATION_RATE);
								if ( VERBOSE ) {
									//System.out.println("flag7: observedBase-" + observedBase + "\t" + "g.base1-" + g.base1 + "\t" + "g.base2-" + g.base2 + "\t" + "likelihood-" + log10(likelihood));
								}
							}
							break;
						case BaseUtils.A:
							if(refPreBase == BaseUtils.C){
								likelihood += observedBase == g.base1 ? pOfBase1 * (1.0-error) : (g.base1 == BaseUtils.G ? pOfBase1 * (error/3.0 + (1.0-CPG_METHYLATION_RATE) * BISULFITE_CONVERSION_RATE) : pOfBase1 * (error/3.0));
								likelihood += observedBase == g.base2 ? pOfBase1 * (1.0-error) : (g.base2 == BaseUtils.G ? pOfBase2 * (error/3.0 + (1.0-CPG_METHYLATION_RATE) * BISULFITE_CONVERSION_RATE) : pOfBase2 * (error/3.0));
								if ( VERBOSE ) {
									//System.out.println("flag8: observedBase-" + observedBase + "\t" + "g.base1-" + g.base1 + "\t" + "g.base2-" + g.base2 + "\t" + "likelihood-" + log10(likelihood));
								}
							}
							else{
								likelihood += observedBase == g.base1 ? pOfBase1 * (1.0-error) : (g.base1 == BaseUtils.G ? pOfBase1 * (error/3.0 + (1.0-CPH_METHYLATION_RATE) * BISULFITE_CONVERSION_RATE) : pOfBase1 * (error/3.0));
								likelihood += observedBase == g.base2 ? pOfBase1 * (1.0-error) : (g.base2 == BaseUtils.G ? pOfBase2 * (error/3.0 + (1.0-CPH_METHYLATION_RATE) * BISULFITE_CONVERSION_RATE) : pOfBase2 * (error/3.0));
								if ( VERBOSE ) {
									//System.out.println("flag9: observedBase-" + observedBase + "\t" + "g.base1-" + g.base1 + "\t" + "g.base2-" + g.base2 + "\t" + "likelihood-" + log10(likelihood));
								}
							}
							
							break;
						default: 
							likelihood += observedBase == g.base1 ? pOfBase1 * (1.0-error) : pOfBase1 * (error/3.0);
							likelihood += observedBase == g.base2 ? pOfBase2 * (1.0-error) : pOfBase2 * (error/3.0);
							if ( VERBOSE ) {
								//System.out.println("flag9: observedBase-" + observedBase + "\t" + "g.base1-" + g.base1 + "\t" + "g.base2-" + g.base2 + "\t" + "likelihood-" + log10(likelihood));
							}
						}
				}
				 double logLikelihood = log10(likelihood);
	             gl.log10Likelihoods[g.ordinal()] += logLikelihood;
	             gl.log10Posteriors[g.ordinal()] += logLikelihood;
			}
			
			if ( VERBOSE ) {
                for ( DiploidGenotype g : DiploidGenotype.values() ) { System.out.printf("%s\t", g); }
                System.out.println();
                for ( DiploidGenotype g : DiploidGenotype.values() ) { System.out.printf("%.2f\t", gl.log10Likelihoods[g.ordinal()]); }
                System.out.println();
                for ( DiploidGenotype g : DiploidGenotype.values() ) { System.out.printf("%.2f\t", gl.getPriors()[g.ordinal()]); }
                System.out.println();
            }

            return gl;
            
		} catch (CloneNotSupportedException e) {
			// TODO Auto-generated catch block
			throw new RuntimeException(e);
		}
		
		
		
	}
	

	protected BisulfiteDiploidSNPGenotypeLikelihoods getCalculateGenotypeLikelihoods(byte observedBase, byte qualityScore, boolean negStrand){

		try {
			BisulfiteDiploidSNPGenotypeLikelihoods gl = (BisulfiteDiploidSNPGenotypeLikelihoods)this.clone();
			gl.setToZeroBs();
			double error = pow(10, -qualityScore/10.0);
			double pOfBase1 = 0.5, pOfBase2 = 1.0 - pOfBase1;
			for ( DiploidGenotype g : DiploidGenotype.values() ) {
				double likelihood = 0.0;
				if(!negStrand){
						switch(observedBase){
							case BaseUtils.C:
								if(!DETERMINED_CYTOSINE_TYPE_NEGSTRAND){
									likelihood += observedBase == g.base1 ? pOfBase1 * (1.0-error) * ( (1.0-DETERMINED_CYTOSINE_TYPE_C_METHY) * (1.0-BISULFITE_CONVERSION_RATE) + DETERMINED_CYTOSINE_TYPE_C_METHY ) : pOfBase1 * (error/3.0) * ( (1.0-DETERMINED_CYTOSINE_TYPE_C_METHY) * (1.0-BISULFITE_CONVERSION_RATE) + DETERMINED_CYTOSINE_TYPE_C_METHY );
									likelihood += observedBase == g.base2 ? pOfBase2 * (1.0-error) *  ( (1.0-DETERMINED_CYTOSINE_TYPE_C_METHY) * (1.0-BISULFITE_CONVERSION_RATE) + DETERMINED_CYTOSINE_TYPE_C_METHY ) : pOfBase2 * (error/3.0) * ( (1.0-DETERMINED_CYTOSINE_TYPE_C_METHY) * (1.0-BISULFITE_CONVERSION_RATE) + DETERMINED_CYTOSINE_TYPE_C_METHY );
								}
								else{
									likelihood += observedBase == g.base1 ? pOfBase1 * (1.0-error) * ( (1.0-CPH_METHYLATION_RATE) * (1.0-BISULFITE_CONVERSION_RATE) + CPH_METHYLATION_RATE ) : pOfBase1 * (error/3.0) * ( (1.0-CPH_METHYLATION_RATE) * (1.0-BISULFITE_CONVERSION_RATE) + CPH_METHYLATION_RATE );
									likelihood += observedBase == g.base2 ? pOfBase2 * (1.0-error) *  ( (1.0-CPH_METHYLATION_RATE) * (1.0-BISULFITE_CONVERSION_RATE) + CPH_METHYLATION_RATE ) : pOfBase2 * (error/3.0) * ( (1.0-CPH_METHYLATION_RATE) * (1.0-BISULFITE_CONVERSION_RATE) + CPH_METHYLATION_RATE );
								}
									
									if ( VERBOSE ) {
									//	System.out.println("flag1: observedBase-" + observedBase + "\t" + "g.base1-" + g.base1 + "\t" + "g.base2-" + g.base2 + "\t" + "likelihood-" + log10(likelihood));
									}
								
								break;
							case BaseUtils.T:
								if(!DETERMINED_CYTOSINE_TYPE_NEGSTRAND){
									likelihood += observedBase == g.base1 ? pOfBase1 * (1.0-error) : (g.base1 == BaseUtils.C ? pOfBase1 * (error/3.0 + (1.0-DETERMINED_CYTOSINE_TYPE_C_METHY) * BISULFITE_CONVERSION_RATE) : pOfBase1 * (error/3.0));
									likelihood += observedBase == g.base2 ? pOfBase1 * (1.0-error) : (g.base2 == BaseUtils.C ? pOfBase2 * (error/3.0 + (1.0-DETERMINED_CYTOSINE_TYPE_C_METHY) * BISULFITE_CONVERSION_RATE) : pOfBase2 * (error/3.0));

								}
								else{
									likelihood += observedBase == g.base1 ? pOfBase1 * (1.0-error) : (g.base1 == BaseUtils.C ? pOfBase1 * (error/3.0 + (1.0-CPH_METHYLATION_RATE) * BISULFITE_CONVERSION_RATE) : pOfBase1 * (error/3.0));
									likelihood += observedBase == g.base2 ? pOfBase1 * (1.0-error) : (g.base2 == BaseUtils.C ? pOfBase2 * (error/3.0 + (1.0-CPH_METHYLATION_RATE) * BISULFITE_CONVERSION_RATE) : pOfBase2 * (error/3.0));

								}
								if ( VERBOSE ) {
									//	System.out.println("flag3: observedBase-" + observedBase + "\t" + "g.base1-" + g.base1 + "\t" + "g.base2-" + g.base2 + "\t" + "likelihood-" + log10(likelihood));
									}
								
								break;
							default: 
								likelihood += observedBase == g.base1 ? pOfBase1 * (1.0-error) : pOfBase1 * (error/3.0);
								likelihood += observedBase == g.base2 ? pOfBase2 * (1.0-error) : pOfBase2 * (error/3.0);
								if ( VERBOSE ) {
									//System.out.println("flag5: observedBase-" + observedBase + "\t" + "g.base1-" + g.base1 + "\t" + "g.base2-" + g.base2 + "\t" + "likelihood-" + log10(likelihood));
								}
									
						}
				}
				else{
						switch(observedBase){
						case BaseUtils.G:
							if(DETERMINED_CYTOSINE_TYPE_NEGSTRAND){
								likelihood += observedBase == g.base1 ? pOfBase1 * (1.0-error) * ( (1.0-DETERMINED_CYTOSINE_TYPE_C_METHY) * (1.0-BISULFITE_CONVERSION_RATE) + DETERMINED_CYTOSINE_TYPE_C_METHY ) : pOfBase1 * (error/3.0) * ( (1.0-DETERMINED_CYTOSINE_TYPE_C_METHY) * (1.0-BISULFITE_CONVERSION_RATE) + DETERMINED_CYTOSINE_TYPE_C_METHY );
								likelihood += observedBase == g.base2 ? pOfBase2 * (1.0-error) *  ( (1.0-DETERMINED_CYTOSINE_TYPE_C_METHY) * (1.0-BISULFITE_CONVERSION_RATE) + DETERMINED_CYTOSINE_TYPE_C_METHY ) : pOfBase2 * (error/3.0) * ( (1.0-DETERMINED_CYTOSINE_TYPE_C_METHY) * (1.0-BISULFITE_CONVERSION_RATE) + DETERMINED_CYTOSINE_TYPE_C_METHY );
								if ( VERBOSE ) {
									//System.out.println("flag6: observedBase-" + observedBase + "\t" + "g.base1-" + g.base1 + "\t" + "g.base2-" + g.base2 + "\t" + "likelihood-" + log10(likelihood));
								}
							}
							else{
								likelihood += observedBase == g.base1 ? pOfBase1 * (1.0-error) * ((1.0-CPH_METHYLATION_RATE) * (1.0-BISULFITE_CONVERSION_RATE) + CPH_METHYLATION_RATE) : pOfBase1 * (error/3.0) * ((1.0-CPH_METHYLATION_RATE) * (1.0-BISULFITE_CONVERSION_RATE) + CPH_METHYLATION_RATE);
								likelihood += observedBase == g.base2 ? pOfBase2 * (1.0-error) * ((1.0-CPH_METHYLATION_RATE) * (1.0-BISULFITE_CONVERSION_RATE) + CPH_METHYLATION_RATE) : pOfBase2 * (error/3.0) * ((1.0-CPH_METHYLATION_RATE) * (1.0-BISULFITE_CONVERSION_RATE) + CPH_METHYLATION_RATE);
								if ( VERBOSE ) {
									//System.out.println("flag7: observedBase-" + observedBase + "\t" + "g.base1-" + g.base1 + "\t" + "g.base2-" + g.base2 + "\t" + "likelihood-" + log10(likelihood));
								}
							}
							break;
						case BaseUtils.A:
							if(DETERMINED_CYTOSINE_TYPE_NEGSTRAND){
								likelihood += observedBase == g.base1 ? pOfBase1 * (1.0-error) : (g.base1 == BaseUtils.G ? pOfBase1 * (error/3.0 + (1.0-DETERMINED_CYTOSINE_TYPE_C_METHY) * BISULFITE_CONVERSION_RATE) : pOfBase1 * (error/3.0));
								likelihood += observedBase == g.base2 ? pOfBase1 * (1.0-error) : (g.base2 == BaseUtils.G ? pOfBase2 * (error/3.0 + (1.0-DETERMINED_CYTOSINE_TYPE_C_METHY) * BISULFITE_CONVERSION_RATE) : pOfBase2 * (error/3.0));
								if ( VERBOSE ) {
									//System.out.println("flag8: observedBase-" + observedBase + "\t" + "g.base1-" + g.base1 + "\t" + "g.base2-" + g.base2 + "\t" + "likelihood-" + log10(likelihood));
								}
							}
							else{
								likelihood += observedBase == g.base1 ? pOfBase1 * (1.0-error) : (g.base1 == BaseUtils.G ? pOfBase1 * (error/3.0 + (1.0-CPH_METHYLATION_RATE) * BISULFITE_CONVERSION_RATE) : pOfBase1 * (error/3.0));
								likelihood += observedBase == g.base2 ? pOfBase1 * (1.0-error) : (g.base2 == BaseUtils.G ? pOfBase2 * (error/3.0 + (1.0-CPH_METHYLATION_RATE) * BISULFITE_CONVERSION_RATE) : pOfBase2 * (error/3.0));
								if ( VERBOSE ) {
									//System.out.println("flag9: observedBase-" + observedBase + "\t" + "g.base1-" + g.base1 + "\t" + "g.base2-" + g.base2 + "\t" + "likelihood-" + log10(likelihood));
								}
							}
							
							break;
						default: 
							likelihood += observedBase == g.base1 ? pOfBase1 * (1.0-error) : pOfBase1 * (error/3.0);
							likelihood += observedBase == g.base2 ? pOfBase2 * (1.0-error) : pOfBase2 * (error/3.0);
							if ( VERBOSE ) {
								//System.out.println("flag9: observedBase-" + observedBase + "\t" + "g.base1-" + g.base1 + "\t" + "g.base2-" + g.base2 + "\t" + "likelihood-" + log10(likelihood));
							}
						}
				}
				 double logLikelihood = log10(likelihood);
	             gl.log10Likelihoods[g.ordinal()] += logLikelihood;
	             gl.log10Posteriors[g.ordinal()] += logLikelihood;
			}
			
			if ( VERBOSE ) {
                for ( DiploidGenotype g : DiploidGenotype.values() ) { System.out.printf("%s\t", g); }
                System.out.println();
                for ( DiploidGenotype g : DiploidGenotype.values() ) { System.out.printf("%.2f\t", gl.log10Likelihoods[g.ordinal()]); }
                System.out.println();
                for ( DiploidGenotype g : DiploidGenotype.values() ) { System.out.printf("%.2f\t", gl.getPriors()[g.ordinal()]); }
                System.out.println();
            }

            return gl;
            
		} catch (CloneNotSupportedException e) {
			// TODO Auto-generated catch block
			throw new RuntimeException(e);
		}
		
		
		
	}
	
    @Override
    public double[] getPriors() {
        return this.priors.getPriors();
    }
    
    public void setPriors(BisulfiteDiploidSNPGenotypePriors priors) {
        this.priors = priors;
        log10Posteriors = new double[DiploidGenotype.values().length];
        for ( DiploidGenotype g : DiploidGenotype.values() ) {
            int i = g.ordinal();
            log10Posteriors[i] = priors.getPriors()[i] + log10Likelihoods[i];
        }
    }
 // likelihoods are all zeros

    protected void setToZeroBs() {
    	//System.err.println("test!!!");
    	for ( DiploidGenotype g : DiploidGenotype.values() ) {
            int i = g.ordinal();
            log10Likelihoods[i] = PROB_OF_REFERENCE_ERROR;
        }

       // if(this.priors.getPriors() == null){
        //	System.err.println("priors null!!!");
       // }
       // System.err.println("log10Likelihoods null!!!");
        log10Posteriors = this.priors.getPriors().clone();     // posteriors are all the priors
    }

    @Override
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
	
}
