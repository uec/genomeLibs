package org.broadinstitute.sting.gatk.bisulfitegenotyper;

import java.util.Map;

import org.apache.log4j.Logger;
import org.broad.tribble.util.variantcontext.Allele;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.StratifiedAlignmentContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.genotyper.BiallelicGenotypeLikelihoods;
import org.broadinstitute.sting.gatk.walkers.genotyper.DiploidSNPGenotypeLikelihoods;
import org.broadinstitute.sting.gatk.walkers.genotyper.DiploidSNPGenotypePriors;
import org.broadinstitute.sting.gatk.walkers.genotyper.GenotypePriors;
import org.broadinstitute.sting.gatk.walkers.genotyper.SNPGenotypeLikelihoodsCalculationModel;
import org.broadinstitute.sting.gatk.walkers.genotyper.UnifiedArgumentCollection;
import org.broadinstitute.sting.gatk.walkers.genotyper.UnifiedGenotyperEngine;
import org.broadinstitute.sting.utils.exceptions.StingException;
import org.broadinstitute.sting.utils.genotype.DiploidGenotype;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;

public class BisulfiteSNPGenotypeLikelihoodsCalculationModel extends
		SNPGenotypeLikelihoodsCalculationModel {

	protected Byte bestAlternateAllele = null;
	protected Byte secondBestAlternateAllele = null;
	
	public BisulfiteSNPGenotypeLikelihoodsCalculationModel(
			UnifiedArgumentCollection UAC, Logger logger) {
		super(UAC, logger);
		// TODO Auto-generated constructor stub
	}
	

	public Allele getLikelihoodsBs(RefMetaDataTracker tracker,
            ReferenceContext ref,
            Map<String, StratifiedAlignmentContext> contexts,
            StratifiedAlignmentContext.StratifiedContextType contextType,
            GenotypePriors priors,
            Map<String, BisulfiteBiallelicGenotypeLikelihoods> GLs,
            Allele alternateAlleleToUse) {
		if ( !(priors instanceof DiploidSNPGenotypePriors) )
            throw new StingException("Only diploid-based SNP priors are supported in the SNP GL model");

        byte refBase = ref.getBase();
        Allele refAllele = Allele.create(refBase, true);
     // if there are no non-ref bases...
        if ( bestAlternateAllele == null ) {
            // if we only want variants, then we don't need to calculate genotype likelihoods
            if ( UAC.OutputMode == UnifiedGenotyperEngine.OUTPUT_MODE.EMIT_VARIANTS_ONLY )
                return refAllele;

            // otherwise, choose any alternate allele (it doesn't really matter)
            bestAlternateAllele = (byte)(refBase != 'A' ? 'A' : 'C');
            secondBestAlternateAllele = (byte)(refBase != 'G' ? 'G' : 'T');
        }

        Allele altAllele = Allele.create(bestAlternateAllele, false);
        Allele secondAltAllele = Allele.create(secondBestAlternateAllele, false);

        for ( Map.Entry<String, StratifiedAlignmentContext> sample : contexts.entrySet() ) {
            ReadBackedPileup pileup = sample.getValue().getContext(contextType).getBasePileup();

            // create the GenotypeLikelihoods object, need change to BisulfiteDiploidSNPGenotypePriors
            BisulfiteDiploidSNPGenotypeLikelihoods GL = new BisulfiteDiploidSNPGenotypeLikelihoods((DiploidSNPGenotypePriors)priors, UAC.PCR_error);
            int nGoodBases = GL.add(pileup, true, true);
            if ( nGoodBases == 0 )
                continue;

            double[] likelihoods = GL.getLikelihoods();

            DiploidGenotype refGenotype = DiploidGenotype.createHomGenotype(refBase);
            DiploidGenotype hetGenotype = DiploidGenotype.createDiploidGenotype(refBase, bestAlternateAllele);
            DiploidGenotype homGenotype = DiploidGenotype.createHomGenotype(bestAlternateAllele);
            DiploidGenotype homNonrefGenotype = DiploidGenotype.createDiploidGenotype(bestAlternateAllele,secondBestAlternateAllele);
            GLs.put(sample.getKey(), new BisulfiteBiallelicGenotypeLikelihoods(sample.getKey(),
                    refAllele,
                    altAllele,
                    secondAltAllele,
                    likelihoods[refGenotype.ordinal()],
                    likelihoods[hetGenotype.ordinal()],
                    likelihoods[homGenotype.ordinal()],
                    likelihoods[homNonrefGenotype.ordinal()],
                    getFilteredDepth(pileup)));
        }

        return refAllele;
	}
	
	
}
