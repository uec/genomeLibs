package org.broadinstitute.sting.gatk.uscec.bisulfitesnpmodel;

import java.util.Map;

import org.apache.log4j.Logger;
import org.broad.tribble.util.variantcontext.Allele;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.StratifiedAlignmentContext;
import org.broadinstitute.sting.gatk.contexts.StratifiedAlignmentContext.StratifiedContextType;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.genotyper.BiallelicGenotypeLikelihoods;
import org.broadinstitute.sting.gatk.walkers.genotyper.GenotypeLikelihoodsCalculationModel;
import org.broadinstitute.sting.gatk.walkers.genotyper.GenotypePriors;
import org.broadinstitute.sting.gatk.walkers.genotyper.UnifiedArgumentCollection;

public class NonRefDependSNPGenotypeLikelihoodsCalculationModel extends
		GenotypeLikelihoodsCalculationModel {

	public enum MethylSNPModel {
        BM,
        GM,
        NM
    }
	
	public NonRefDependSNPGenotypeLikelihoodsCalculationModel(
			UnifiedArgumentCollection UAC, Logger logger) {
		super(UAC, logger);
		// TODO Auto-generated constructor stub
	}

	@Override
	public Allele getLikelihoods(RefMetaDataTracker tracker,
			ReferenceContext ref,
			Map<String, StratifiedAlignmentContext> contexts,
			StratifiedContextType contextType, GenotypePriors priors,
			Map<String, BiallelicGenotypeLikelihoods> GLs,
			Allele alternateAlleleToUse) {
		// TODO Auto-generated method stub
		return null;
	}
	
	public Allele getBsLikelihoods(RefMetaDataTracker tracker,
			ReferenceContext ref,
			Map<String, StratifiedAlignmentContext> contexts,
			StratifiedContextType contextType,
			Map<String, BisulfiteBiallelicGenotypeLikelihoods> GLs,
			Allele alternateAlleleToUse) {
		// TODO Auto-generated method stub
		return null;
	}
	
	public void initialize(CytosineTypeStatus cts, BisulfiteArgumentCollection BAC, boolean autoEstimateC, boolean secondIteration){

	}

}
