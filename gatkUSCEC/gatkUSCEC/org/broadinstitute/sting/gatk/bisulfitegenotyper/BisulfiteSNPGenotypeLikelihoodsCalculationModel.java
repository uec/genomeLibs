package org.broadinstitute.sting.gatk.bisulfitegenotyper;

import java.util.Map;

import net.sf.samtools.SAMRecord;

import org.apache.log4j.Logger;
import org.broad.tribble.util.variantcontext.Allele;
import org.broad.tribble.util.variantcontext.VariantContext;
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
import org.broadinstitute.sting.gatk.walkers.genotyper.GenotypeLikelihoodsCalculationModel.GENOTYPING_MODE;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.exceptions.StingException;
import org.broadinstitute.sting.utils.genotype.DiploidGenotype;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;

public class BisulfiteSNPGenotypeLikelihoodsCalculationModel extends
		SNPGenotypeLikelihoodsCalculationModel {

	protected Byte bestAlternateAllele = null;
	protected Byte secondBestAlternateAllele = null;
	protected final boolean useAlleleFromVCF;
	protected int testLoc = 7152084;
	
	public BisulfiteSNPGenotypeLikelihoodsCalculationModel(
			UnifiedArgumentCollection UAC, Logger logger) {
		super(UAC, logger);
		useAlleleFromVCF = UAC.GenotypingMode == GENOTYPING_MODE.GENOTYPE_GIVEN_ALLELES;
		// TODO Auto-generated constructor stub
	}
	

	public Allele getLikelihoodsBs(RefMetaDataTracker tracker,
            ReferenceContext ref,
            Map<String, StratifiedAlignmentContext> contexts,
            StratifiedAlignmentContext.StratifiedContextType contextType,
            GenotypePriors priors,
            Map<String, BisulfiteBiallelicGenotypeLikelihoods> GLs,
            Allele alternateAlleleToUse) {
		if ( !(priors instanceof BisulfiteDiploidSNPGenotypePriors) )
            throw new StingException("Only Bisulfite diploid-based SNP priors are supported in the SNP GL model");

        byte refBase = ref.getBase();
        Allele refAllele = Allele.create(refBase, true);
        

        // find the alternate allele with the largest sum of quality scores
        if ( alternateAlleleToUse != null ) {
            bestAlternateAllele = alternateAlleleToUse.getBases()[0];
            secondBestAlternateAllele = alternateAlleleToUse.getBases()[0];
        } else if ( useAlleleFromVCF ) {
            final VariantContext vcInput = tracker.getVariantContext(ref, "alleles", null, ref.getLocus(), true);
            if ( vcInput == null )
                return null;
            if ( !vcInput.isSNP() ) {
                logger.info("Record at position " + ref.getLocus() + " is not a SNP; skipping...");
                return null;
            }
            if ( !vcInput.isBiallelic() ) {
                logger.info("Record at position " + ref.getLocus() + " is not bi-allelic; choosing the first allele...");
                //return null;
            }
            bestAlternateAllele = vcInput.getAlternateAllele(0).getBases()[0];
            secondBestAlternateAllele = vcInput.getAlternateAllele(0).getBases()[0];
        } else {
            initializeBestAlternateAllele(refBase, contexts);

        }
        
      
        //System.err.println("ok");
     // if there are no non-ref bases...
        if ( bestAlternateAllele == null) {
            // if we only want variants, then we don't need to calculate genotype likelihoods
            if ( UAC.OutputMode == UnifiedGenotyperEngine.OUTPUT_MODE.EMIT_VARIANTS_ONLY )
                return refAllele;

            // otherwise, choose any alternate allele (it doesn't really matter)
            bestAlternateAllele = (byte)(refBase != BaseUtils.A ? BaseUtils.A : BaseUtils.C);
            secondBestAlternateAllele = (byte)(refBase != BaseUtils.G ? BaseUtils.G : BaseUtils.T);
        }
        
        

        Allele altAllele = Allele.create(bestAlternateAllele, false);
        Allele secondAltAllele = altAllele;
        if(secondBestAlternateAllele != null){
        	secondAltAllele = Allele.create(secondBestAlternateAllele, false);
        }
        

        for ( Map.Entry<String, StratifiedAlignmentContext> sample : contexts.entrySet() ) {
            ReadBackedPileup pileup = sample.getValue().getContext(contextType).getBasePileup();

            // do not use this prior, this prior is flat prior intiated in genotypeEngine, so we actually do not transfer this priors...
            BisulfiteDiploidSNPGenotypeLikelihoods GL = new BisulfiteDiploidSNPGenotypeLikelihoods(tracker, ref, (BisulfiteDiploidSNPGenotypePriors)priors, UAC.PCR_error);
            int nGoodBases = GL.add(pileup, true, true);
            if ( nGoodBases == 0 )
                continue;

            double[] likelihoods = GL.getLikelihoods();
            double[] posterir = GL.getPosteriors();
            double[] prio = GL.getPriors();
            

            DiploidGenotype refGenotype = DiploidGenotype.createHomGenotype(refBase);
            DiploidGenotype hetGenotype = DiploidGenotype.createDiploidGenotype(refBase, bestAlternateAllele);
            DiploidGenotype homGenotype = DiploidGenotype.createHomGenotype(bestAlternateAllele);
            if((pileup.getLocation().getStart()) == testLoc){
            	System.out.println("sample: " + sample.getKey());
            	System.out.println("sample location: " + pileup.getPileupString((char)refBase));
            	System.out.println("sample: " + sample.getValue().getLocation().getStart());
            	System.out.println("refBase: " + refBase + " bestAlternateAllele: " + bestAlternateAllele);
            	System.out.println("nGoodBases " + nGoodBases);
            	System.out.println("refGenotype " + likelihoods[refGenotype.ordinal()] + "\t" + prio[refGenotype.ordinal()] + "\t" + posterir[refGenotype.ordinal()]);
            	System.out.println("hetGenotype " + likelihoods[hetGenotype.ordinal()] + "\t" + prio[hetGenotype.ordinal()] + "\t" + posterir[hetGenotype.ordinal()]);
            	System.out.println("homGenotype " + likelihoods[homGenotype.ordinal()] + "\t" + prio[homGenotype.ordinal()] + "\t" + posterir[homGenotype.ordinal()]);
            }
            if(secondBestAlternateAllele != null){
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
            else{
            	GLs.put(sample.getKey(), new BisulfiteBiallelicGenotypeLikelihoods(sample.getKey(),
                        refAllele,
                        altAllele,
                        likelihoods[refGenotype.ordinal()],
                        likelihoods[hetGenotype.ordinal()],
                        likelihoods[homGenotype.ordinal()],
                        getFilteredDepth(pileup)));
            }
            
        }

        return refAllele;
	}
	
	@Override
	protected void initializeBestAlternateAllele(byte ref, Map<String, StratifiedAlignmentContext> contexts) {
        int[] qualCounts = new int[4];
        int location = 0;
        for ( Map.Entry<String, StratifiedAlignmentContext> sample : contexts.entrySet() ) {
            // calculate the sum of quality scores for each base
            ReadBackedPileup pileup = sample.getValue().getContext(StratifiedAlignmentContext.StratifiedContextType.COMPLETE).getBasePileup();
            for ( PileupElement p : pileup ) {
                // ignore deletions and filtered bases
            	SAMRecord samRecord = p.getRead();
            	
            	boolean negStrand = samRecord.getReadNegativeStrandFlag();
				int alignmentS = samRecord.getAlignmentStart();
				int	onRefCoord = (negStrand) ? samRecord.getUnclippedEnd() : alignmentS;

				if((pileup.getLocation().getStart()) == testLoc){
					System.out.println("before filter:\t" + onRefCoord + "\t" + p.getOffset() + "\t" + negStrand + "\t" + pileup.getLocation().getStart() + "\t" + (char)p.getBase());
					System.out.println("deletion: " + p.isDeletion());
					System.out.println("GATKSAMRecord: " + (p.getRead() instanceof GATKSAMRecord));
					System.out.println("isGoodBase: " + ((GATKSAMRecord)p.getRead()).isGoodBase(p.getOffset()));
		                     
				}
                if ( p.isDeletion() ||
                     (p.getRead() instanceof GATKSAMRecord && !((GATKSAMRecord)p.getRead()).isGoodBase(p.getOffset())) )
                    continue;
                //if(p.getBase() != ref){
                	if((p.getRead().getReadNegativeStrandFlag() && p.getBase() == BaseUtils.A) || (!p.getRead().getReadNegativeStrandFlag() && p.getBase() == BaseUtils.T))
                		continue;
                //}
                	
                	
					if((pileup.getLocation().getStart()) == testLoc){
						System.out.println(onRefCoord + "\t" + p.getOffset() + "\t" + negStrand + "\t" + pileup.getLocation().getStart() + "\t" + (char)p.getBase());
						location = pileup.getLocation().getStart();
					}
                int index = BaseUtils.simpleBaseToBaseIndex(p.getBase());
                if ( index >= 0 )
                    qualCounts[index] += p.getQual();
            }
        }

        // set the non-ref base with maximum quality score sum
        int maxCount = 0;
        int secondMaxCount = 0;
        bestAlternateAllele = null;
        secondBestAlternateAllele = null;
        for ( byte altAllele : BaseUtils.BASES ) {
            if ( altAllele == ref )
                continue;
            //if( altAllele == 'T')
            	//continue;
            int index = BaseUtils.simpleBaseToBaseIndex(altAllele);
            if ( qualCounts[index] > maxCount ) {
            	secondMaxCount = maxCount;
            	maxCount = qualCounts[index];
            	if(bestAlternateAllele != null){
            		secondBestAlternateAllele = bestAlternateAllele;
            	}
                bestAlternateAllele = altAllele;
            }
            else if (qualCounts[index] > secondMaxCount && qualCounts[index] < maxCount){
            	secondMaxCount = qualCounts[index];
            	secondBestAlternateAllele = altAllele;
            }
        }
        
        if(location == testLoc){
        	System.out.println("bestAlternateAllele: " + bestAlternateAllele + "\t" + maxCount);
        	if(secondBestAlternateAllele != null){
        		System.out.println("secondBestAlternateAllele: " + "\t" + secondBestAlternateAllele + "\t" + secondMaxCount);
        	}
        }

    }
	
	
	
	
}
