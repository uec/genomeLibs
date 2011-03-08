package org.broadinstitute.sting.gatk.bisulfitegenotyper;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;

import org.broad.tribble.vcf.VCFConstants;
import org.broad.tribble.vcf.VCFFilterHeaderLine;
import org.broad.tribble.vcf.VCFHeader;
import org.broad.tribble.vcf.VCFHeaderLine;
import org.broad.tribble.vcf.VCFHeaderLineType;
import org.broad.tribble.vcf.VCFInfoHeaderLine;
import org.broad.tribble.vcf.VCFWriter;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.ArgumentCollection;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.DownsampleType;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.datasources.rmd.ReferenceOrderedDataSource;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.utils.helpers.DbSNPHelper;
import org.broadinstitute.sting.gatk.walkers.BAQMode;
import org.broadinstitute.sting.gatk.walkers.By;
import org.broadinstitute.sting.gatk.walkers.DataSource;
import org.broadinstitute.sting.gatk.walkers.Downsample;
import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.gatk.walkers.Reference;
import org.broadinstitute.sting.gatk.walkers.TreeReducible;
import org.broadinstitute.sting.gatk.walkers.Window;
import org.broadinstitute.sting.gatk.walkers.annotator.VariantAnnotatorEngine;
import org.broadinstitute.sting.gatk.walkers.genotyper.GenotypeLikelihoodsCalculationModel;
import org.broadinstitute.sting.gatk.walkers.genotyper.UnifiedArgumentCollection;
import org.broadinstitute.sting.gatk.walkers.genotyper.UnifiedGenotyper;
import org.broadinstitute.sting.gatk.walkers.genotyper.UnifiedGenotyperEngine;
import org.broadinstitute.sting.gatk.walkers.genotyper.VariantCallContext;
import org.broadinstitute.sting.utils.SampleUtils;
import org.broadinstitute.sting.utils.baq.BAQ;
import org.broadinstitute.sting.utils.vcf.VCFUtils;



/**
* A variant caller which unifies the approaches of several disparate callers.  Works for single-sample and
* multi-sample data.  The user can choose from several different incorporated calculation models.
*/
//todo -- change when UG is generalized to do BAQ as necessary
//@BAQMode(QualityMode = BAQ.QualityMode.DONT_MODIFY, ApplicationTime = BAQ.ApplicationTime.HANDLED_IN_WALKER)
@BAQMode(QualityMode = BAQ.QualityMode.OVERWRITE_QUALS, ApplicationTime = BAQ.ApplicationTime.ON_INPUT)
@Reference(window=@Window(start=-200,stop=200))
@By(DataSource.REFERENCE)
@Downsample(by=DownsampleType.BY_SAMPLE, toCoverage=250)
public class BisulfiteGenotyper extends LocusWalker<VariantCallContext, BisulfiteGenotyper.BGStatistics> implements TreeReducible<BisulfiteGenotyper.BGStatistics> {

	
	
	@ArgumentCollection protected UnifiedArgumentCollection UAC = new UnifiedArgumentCollection();

    // control the output
    @Output(doc="File to which variants should be written",required=true)
    protected VCFWriter writer = null;

    @Argument(fullName = "verbose_mode", shortName = "verbose", doc = "File to print all of the annotated and detailed debugging output", required = false)
    protected PrintStream verboseWriter = null;

    @Argument(fullName = "metrics_file", shortName = "metrics", doc = "File to print any relevant callability metrics output", required = false)
    protected PrintStream metricsWriter = null;

    @Argument(fullName="annotation", shortName="A", doc="One or more specific annotations to apply to variant calls", required=false)
    protected List<String> annotationsToUse = new ArrayList<String>();

    @Argument(fullName="group", shortName="G", doc="One or more classes/groups of annotations to apply to variant calls", required=false)
    protected String[] annotationClassesToUse = { "Standard" };

    // the calculation arguments
    private BisulfiteGenotyperEngine BG_engine = null;

    // the annotation engine
    protected VariantAnnotatorEngine annotationEngine;

    // enable deletions in the pileup
    public boolean includeReadsWithDeletionAtLoci() { return true; }

    // enable extended events for indels
    public boolean generateExtendedEvents() { return UAC.GLmodel == GenotypeLikelihoodsCalculationModel.Model.DINDEL; }
	
	/**
     * Inner class for collecting output statistics from the UG
     */
    public static class BGStatistics {
        /** The total number of passes examined -- i.e., the number of map calls */
        long nBasesVisited = 0;

        /** The number of bases that were potentially callable -- i.e., those not at excessive coverage or masked with N */
        long nBasesCallable = 0;

        /** The number of bases called confidently (according to user threshold), either ref or other */
        long nBasesCalledConfidently = 0;

        /** The number of bases for which calls were emitted */
        long nCallsMade = 0;

        /** The total number of extended events encountered */
        long nExtendedEvents = 0;

        double percentCallableOfAll()    { return (100.0 * nBasesCallable) / (nBasesVisited-nExtendedEvents); }
        double percentCalledOfAll()      { return (100.0 * nBasesCalledConfidently) / (nBasesVisited-nExtendedEvents); }
        double percentCalledOfCallable() { return (100.0 * nBasesCalledConfidently) / (nBasesCallable); }
    }
    


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
	
	

   

    protected Set<VCFHeaderLine> getHeaderInfo() {
        Set<VCFHeaderLine> headerInfo = new HashSet<VCFHeaderLine>();

        // all annotation fields from VariantAnnotatorEngine
        headerInfo.addAll(annotationEngine.getVCFAnnotationDescriptions());

        // annotation (INFO) fields from UnifiedGenotyper
        if ( !UAC.NO_SLOD )
            headerInfo.add(new VCFInfoHeaderLine(VCFConstants.STRAND_BIAS_KEY, 1, VCFHeaderLineType.Float, "Strand Bias"));
        headerInfo.add(new VCFInfoHeaderLine(VCFConstants.DOWNSAMPLED_KEY, 0, VCFHeaderLineType.Flag, "Were any of the samples downsampled?"));

        // also, check to see whether comp rods were included
        List<ReferenceOrderedDataSource> dataSources = getToolkit().getRodDataSources();
        for ( ReferenceOrderedDataSource source : dataSources ) {
            if ( source.getName().equals(DbSNPHelper.STANDARD_DBSNP_TRACK_NAME) ) {
                headerInfo.add(new VCFInfoHeaderLine(VCFConstants.DBSNP_KEY, 0, VCFHeaderLineType.Flag, "dbSNP Membership"));
            }
            else if ( source.getName().startsWith(VariantAnnotatorEngine.dbPrefix) ) {
                String name = source.getName().substring(VariantAnnotatorEngine.dbPrefix.length());
                headerInfo.add(new VCFInfoHeaderLine(name, 0, VCFHeaderLineType.Flag, name + " Membership"));
            }
        }

        // FORMAT and INFO fields
        headerInfo.addAll(VCFUtils.getSupportedHeaderStrings(VCFConstants.PHRED_GENOTYPE_LIKELIHOODS_KEY));

        // FILTER fields
        if ( UAC.STANDARD_CONFIDENCE_FOR_EMITTING < UAC.STANDARD_CONFIDENCE_FOR_CALLING )
            headerInfo.add(new VCFFilterHeaderLine(UnifiedGenotyperEngine.LOW_QUAL_FILTER_NAME, "Low quality"));

        return headerInfo;
    }

   

    public BGStatistics reduceInit() { return new BGStatistics(); }

    public BGStatistics treeReduce(BGStatistics lhs, BGStatistics rhs) {
        lhs.nBasesCallable += rhs.nBasesCallable;
        lhs.nBasesCalledConfidently += rhs.nBasesCalledConfidently;
        lhs.nBasesVisited += rhs.nBasesVisited;
        lhs.nCallsMade += rhs.nCallsMade;
        return lhs;
    }

    public BGStatistics reduce(VariantCallContext value, BGStatistics sum) {
        // we get a point for reaching reduce
        sum.nBasesVisited++;

        // can't call the locus because of no coverage
        if ( value == null )
            return sum;

        // A call was attempted -- the base was potentially callable
        sum.nBasesCallable++;

        // the base was confidently callable
        sum.nBasesCalledConfidently += value.confidentlyCalled ? 1 : 0;

        // can't make a confident variant call here
        if ( value.vc == null ){
			//System.out.println("No vc ");
			return sum;
        }
            

        try {
            // we are actually making a call
			//System.out.println("ok "+ value.refBase);
            sum.nCallsMade++;
            writer.add(value.vc, value.refBase);
        } catch (IllegalArgumentException e) {
            throw new IllegalArgumentException(e.getMessage() + "; this is often caused by using the --assume_single_sample_reads argument with the wrong sample name");
        }

        return sum;
    }

    public void onTraversalDone(BGStatistics sum) {
        logger.info(String.format("Visited bases                                %d", sum.nBasesVisited));
        logger.info(String.format("Callable bases                               %d", sum.nBasesCallable));
        logger.info(String.format("Confidently called bases                     %d", sum.nBasesCalledConfidently));
        logger.info(String.format("%% callable bases of all loci                 %3.3f", sum.percentCallableOfAll()));
        logger.info(String.format("%% confidently called bases of all loci       %3.3f", sum.percentCalledOfAll()));
        logger.info(String.format("%% confidently called bases of callable loci  %3.3f", sum.percentCalledOfCallable()));
        logger.info(String.format("Actual calls made                            %d", sum.nCallsMade));

        if ( metricsWriter != null ) {
            metricsWriter.println(String.format("Visited bases                                %d", sum.nBasesVisited));
            metricsWriter.println(String.format("Callable bases                               %d", sum.nBasesCallable));
            metricsWriter.println(String.format("Confidently called bases                     %d", sum.nBasesCalledConfidently));
            metricsWriter.println(String.format("%% callable bases of all loci                 %3.3f", sum.percentCallableOfAll()));
            metricsWriter.println(String.format("%% confidently called bases of all loci       %3.3f", sum.percentCalledOfAll()));
            metricsWriter.println(String.format("%% confidently called bases of callable loci  %3.3f", sum.percentCalledOfCallable()));
            metricsWriter.println(String.format("Actual calls made                            %d", sum.nCallsMade));
        }
    }
}