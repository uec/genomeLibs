package org.broadinstitute.sting.gatk.bisulfitegenotyper;

import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.log4j.Logger;
import org.broad.tribble.util.variantcontext.Allele;
import org.broad.tribble.util.variantcontext.Genotype;
import org.broad.tribble.util.variantcontext.GenotypeLikelihoods;
import org.broad.tribble.util.variantcontext.VariantContext;
import org.broad.tribble.vcf.VCFConstants;
import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.StratifiedAlignmentContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.annotator.VariantAnnotatorEngine;
import org.broadinstitute.sting.gatk.walkers.genotyper.BiallelicGenotypeLikelihoods;
import org.broadinstitute.sting.gatk.walkers.genotyper.DindelGenotypeLikelihoodsCalculationModel;
import org.broadinstitute.sting.gatk.walkers.genotyper.DiploidIndelGenotypePriors;
import org.broadinstitute.sting.gatk.walkers.genotyper.DiploidSNPGenotypePriors;
import org.broadinstitute.sting.gatk.walkers.genotyper.GenotypeLikelihoodsCalculationModel;
import org.broadinstitute.sting.gatk.walkers.genotyper.GenotypePriors;
import org.broadinstitute.sting.gatk.walkers.genotyper.SNPGenotypeLikelihoodsCalculationModel;
import org.broadinstitute.sting.gatk.walkers.genotyper.UnifiedArgumentCollection;
import org.broadinstitute.sting.gatk.walkers.genotyper.UnifiedGenotyperEngine;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.fasta.CachingIndexedFastaSequenceFile;

import org.broadinstitute.sting.gatk.bisulfitegenotyper.BisulfiteSNPGenotypeLikelihoodsCalculationModel;

public class BisulfiteGenotyperEngine extends UnifiedGenotyperEngine {

	protected ThreadLocal<BisulfiteSNPGenotypeLikelihoodsCalculationModel> glcm = new ThreadLocal<BisulfiteSNPGenotypeLikelihoodsCalculationModel>();
	
	public BisulfiteGenotyperEngine(GenomeAnalysisEngine toolkit,
			UnifiedArgumentCollection UAC) {
		super(toolkit, UAC);
		// TODO Auto-generated constructor stub
	}

	public BisulfiteGenotyperEngine(GenomeAnalysisEngine toolkit,
			UnifiedArgumentCollection UAC, Logger logger,
			PrintStream verboseWriter, VariantAnnotatorEngine engine,
			Set<String> samples) {
		super(toolkit, UAC, logger, verboseWriter, engine, samples);
		// TODO Auto-generated constructor stub
	}
	
	@Override
	protected void initialize(GenomeAnalysisEngine toolkit, UnifiedArgumentCollection UAC, Logger logger, PrintStream verboseWriter, VariantAnnotatorEngine engine, int numSamples) {
        // note that, because we cap the base quality by the mapping quality, minMQ cannot be less than minBQ
        this.UAC = UAC.clone();
        this.UAC.MIN_MAPPING_QUALTY_SCORE = Math.max(UAC.MIN_MAPPING_QUALTY_SCORE, UAC.MIN_BASE_QUALTY_SCORE);

        this.logger = logger;
        this.verboseWriter = verboseWriter;
        this.annotationEngine = engine;

        N = 2 * numSamples;
        log10AlleleFrequencyPriors = new double[N+1];
        computeAlleleFrequencyPriors(N);
        genotypePriors = this.createGenotypePriors(UAC);

        filter.add(LOW_QUAL_FILTER_NAME);

        try {
            referenceReader = new CachingIndexedFastaSequenceFile(toolkit.getArguments().referenceFile);
        }
        catch(FileNotFoundException ex) {
            throw new UserException.CouldNotReadInputFile(toolkit.getArguments().referenceFile,ex);
        }
    }
	
	protected static GenotypePriors createGenotypePriors(UnifiedArgumentCollection UAC) {
        GenotypePriors priors;
        switch ( UAC.GLmodel ) {
            case SNP:
                // use flat priors for GLs
                priors = new BisulfiteDiploidSNPGenotypePriors();
                break;
            case DINDEL:
                // create flat priors for Indels, actual priors will depend on event length to be genotyped
                priors = new DiploidIndelGenotypePriors();
                break;
            default: throw new IllegalArgumentException("Unexpected GenotypeCalculationModel " + UAC.GLmodel);
        }

        return priors;
    }
	
	@Override
	/**
     * Compute GLs at a given locus.
     *
     * @param tracker    the meta data tracker
     * @param refContext the reference base
     * @param rawContext contextual information around the locus
     * @param alternateAlleleToUse the alternate allele to use, null if not set
     * @return the VariantContext object
     */
    public VariantContext calculateLikelihoods(RefMetaDataTracker tracker, ReferenceContext refContext, AlignmentContext rawContext, Allele alternateAlleleToUse) {
        Map<String, StratifiedAlignmentContext> stratifiedContexts = getFilteredAndStratifiedContexts(UAC, refContext, rawContext);
        VariantContext vc = calculateLikelihoods(tracker, refContext, stratifiedContexts, StratifiedAlignmentContext.StratifiedContextType.COMPLETE, alternateAlleleToUse);
        return GLsToPLs(vc);
    }

	@Override
    protected VariantContext calculateLikelihoods(RefMetaDataTracker tracker, ReferenceContext refContext, Map<String, StratifiedAlignmentContext> stratifiedContexts, StratifiedAlignmentContext.StratifiedContextType type, Allele alternateAlleleToUse) {
		if ( stratifiedContexts == null ){
			//System.out.println("no stratifiedContexts now");
			 return null;
		}

        // initialize the data for this thread if that hasn't been done yet
        if ( glcm.get() == null ) {
            glcm.set(this.getGenotypeLikelihoodsCalculationObject(logger, UAC));
        }

        Map<String, BisulfiteBiallelicGenotypeLikelihoods> GLs = new HashMap<String, BisulfiteBiallelicGenotypeLikelihoods>();

        Allele refAllele = glcm.get().getLikelihoodsBs(tracker, refContext, stratifiedContexts, type, genotypePriors, GLs, alternateAlleleToUse);

        if (refAllele != null){
        	//System.out.println("there is refAllele now");
        	return createVariantContextFromLikelihoodsBs(refContext, refAllele, GLs);
        }   
        else{
        	//System.out.println("no refAllele now");
			 return null;
        }
    }
	

	protected static BisulfiteSNPGenotypeLikelihoodsCalculationModel getGenotypeLikelihoodsCalculationObject(Logger logger, UnifiedArgumentCollection UAC) {
		BisulfiteSNPGenotypeLikelihoodsCalculationModel glcm;
        switch ( UAC.GLmodel ) {
            case SNP:
                glcm = new BisulfiteSNPGenotypeLikelihoodsCalculationModel(UAC, logger);
                break;
           // case DINDEL:
            //    glcm = new DindelGenotypeLikelihoodsCalculationModel(UAC, logger);
            //    break;
            default: throw new IllegalArgumentException("Unexpected GenotypeCalculationModel " + UAC.GLmodel);
        }

        return glcm;
    }

	protected VariantContext createVariantContextFromLikelihoodsBs(ReferenceContext refContext, Allele refAllele, Map<String, BisulfiteBiallelicGenotypeLikelihoods> GLs) {
        // no-call everyone for now
        List<Allele> noCall = new ArrayList<Allele>();
        noCall.add(Allele.NO_CALL);

        Set<Allele> alleles = new HashSet<Allele>();
        alleles.add(refAllele);
        boolean addedAltAllele = false;

        HashMap<String, Genotype> genotypes = new HashMap<String, Genotype>();
        for ( BisulfiteBiallelicGenotypeLikelihoods GL : GLs.values() ) {
            if ( !addedAltAllele ) {
                addedAltAllele = true;
                alleles.add(GL.getAlleleA());
                alleles.add(GL.getAlleleB());
            }

            HashMap<String, Object> attributes = new HashMap<String, Object>();
            GenotypeLikelihoods likelihoods = new GenotypeLikelihoods(GL.getLikelihoods());
            attributes.put(VCFConstants.DEPTH_KEY, GL.getDepth());
            attributes.put(VCFConstants.GENOTYPE_LIKELIHOODS_KEY, likelihoods);

            genotypes.put(GL.getSample(), new Genotype(GL.getSample(), noCall, Genotype.NO_NEG_LOG_10PERROR, null, attributes, false));
        }

        GenomeLoc loc = refContext.getLocus();
        int endLoc = calculateEndPos(alleles, refAllele, loc);

        return new VariantContext("BG_call",
                loc.getContig(),
                loc.getStart(),
                endLoc,
                alleles,
                genotypes,
                VariantContext.NO_NEG_LOG_10PERROR,
                null,
                null);
    }
	
	
	
}
