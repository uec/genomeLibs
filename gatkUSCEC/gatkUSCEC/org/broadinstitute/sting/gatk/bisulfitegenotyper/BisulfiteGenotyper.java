package org.broadinstitute.sting.gatk.bisulfitegenotyper;

import java.util.Arrays;
import java.util.Set;
import java.util.TreeSet;

import org.broad.tribble.vcf.VCFHeader;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.annotator.VariantAnnotatorEngine;
import org.broadinstitute.sting.gatk.walkers.genotyper.UnifiedGenotyper;
import org.broadinstitute.sting.gatk.walkers.genotyper.UnifiedGenotyperEngine;
import org.broadinstitute.sting.gatk.walkers.genotyper.VariantCallContext;
import org.broadinstitute.sting.utils.SampleUtils;

public class BisulfiteGenotyper extends UnifiedGenotyper {

	private UnifiedGenotyperEngine BG_engine = null;
	
	public BisulfiteGenotyper() {
		// TODO Auto-generated constructor stub
	}

	@Override
	public void initialize() {
        // get all of the unique sample names
        // if we're supposed to assume a single sample, do so
        Set<String> samples = new TreeSet<String>();
        if ( UAC.ASSUME_SINGLE_SAMPLE != null )
            samples.add(UAC.ASSUME_SINGLE_SAMPLE);
        else
            samples = SampleUtils.getSAMFileSamples(getToolkit().getSAMFileHeader());

        // initialize the verbose writer
        if ( verboseWriter != null )
            verboseWriter.println("AFINFO\tLOC\tREF\tALT\tMAF\tF\tAFprior\tAFposterior\tNormalizedPosterior");

        annotationEngine = new VariantAnnotatorEngine(getToolkit(), Arrays.asList(annotationClassesToUse), annotationsToUse);
        BG_engine = new BisulfiteGenotyperEngine(getToolkit(), UAC, logger, verboseWriter, annotationEngine, samples);

        // initialize the header
        writer.writeHeader(new VCFHeader(getHeaderInfo(), samples)) ;
    }
	
	@Override
	/**
     * Compute at a given locus.
     *
     * @param tracker the meta data tracker
     * @param refContext the reference base
     * @param rawContext contextual information around the locus
     * @return the VariantCallContext object
     */
    public VariantCallContext map(RefMetaDataTracker tracker, ReferenceContext refContext, AlignmentContext rawContext) {
        return BG_engine.calculateLikelihoodsAndGenotypes(tracker, refContext, rawContext);
    }
}
