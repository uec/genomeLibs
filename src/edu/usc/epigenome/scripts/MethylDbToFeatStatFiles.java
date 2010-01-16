package edu.usc.epigenome.scripts;

import java.io.File;
import java.io.FileOutputStream;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;

import org.biojava.bio.program.gff.GFFRecord;
import org.biojava.bio.program.gff.SimpleGFFRecord;
import org.biojava.bio.seq.StrandedFeature;

import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;
import org.usckeck.genome.ChromFeatures;
import org.usckeck.genome.GFFUtils;

import com.googlecode.charts4j.Color;


import edu.usc.epigenome.genomeLibs.Charts4jUtils;
import edu.usc.epigenome.genomeLibs.ListUtils;
import edu.usc.epigenome.genomeLibs.MatUtils;
import edu.usc.epigenome.genomeLibs.MiscUtils;
import edu.usc.epigenome.genomeLibs.FeatAligners.FeatAligner;
import edu.usc.epigenome.genomeLibs.FeatAligners.FeatAlignerEachfeat;
import edu.usc.epigenome.genomeLibs.GenomicRange.GenomicRangeWithRefpoint;
import edu.usc.epigenome.genomeLibs.MethylDb.Cpg;
import edu.usc.epigenome.genomeLibs.MethylDb.CpgIteratorMultisample;
import edu.usc.epigenome.genomeLibs.MethylDb.MethylDbQuerier;
import edu.usc.epigenome.genomeLibs.MethylDb.MethylDbUtils;
import edu.usc.epigenome.genomeLibs.MethylDb.CpgSummarizers.CpgCoverageSummarizer;
import edu.usc.epigenome.genomeLibs.MethylDb.CpgSummarizers.CpgMethLevelSummarizer;
import edu.usc.epigenome.genomeLibs.MethylDb.CpgSummarizers.CpgSummarizer;


public class MethylDbToFeatStatFiles {

	private static final String C_USAGE = "Use: MethylDbToFeatStatFiles -featShoreStart 500 -featShoreEnd 5000 -onechrom -minCTreads 0 -maxOppStrandAfrac 0.1 -maxFeatSize 10 -centeredSize 0 " +
	" sample1_tablePrefix sample2_tablePrefix ... , feats1.gtf feats2.gtf ...";
	
	@Option(name="-featFilter",usage="We will take the intersection with this feature. Must be a featType in the features table")
	protected List<String> featFilters = new ArrayList<String>(5);
	@Option(name="-maxFeatSize",usage="maximum size of features to include (default Inf)")
    protected int maxFeatSize = Integer.MAX_VALUE;
    @Option(name="-centeredSize",usage="Instead of using the size of the feature, use a set size centered on the feature (default 0=off)")
    protected int centeredSize = 0;
	@Option(name="-onechrom",usage="If set, we only look at chromosome 1 (use for huge feature sets)")
	protected boolean onechrom = false;
    @Option(name="-minCTreads",usage="Minimum number of C or T reads to count as a methylation value")
    protected int minCTreads = 0;
    @Option(name="-featShoreStart",usage="If set, (only) the shore of the feature is used, from featShoreStart bases away to out to featShoreEnd bp")
    protected int featShoreStart = 0;
    @Option(name="-featShoreEnd",usage="If set, (only) the shore of the feature is used, from featShoreStart bases away to out to featShoreEnd bp")
    protected int featShoreEnd = 0;
    @Option(name="-maxOppStrandAfrac",usage="As on the opposite strand are evidence for mutation or SNP. " +
    		"This sets a maximum number of observed As on the opposite strand (default 0.1)")
    protected double maxOppStrandAfrac = 0.1;
	// receives other command line parameters than options
	@Argument
	private List<String> arguments = new ArrayList<String>();

	protected CpgSummarizer[] methSummarizers = new CpgSummarizer[10]; 
	protected CpgSummarizer[] coverageSummarizers = new CpgSummarizer[10];


	/**
	 * @param args
	 */
	public static void main(String[] args)
	throws Exception	
	{
		Logger.getLogger(Logger.GLOBAL_LOGGER_NAME).setLevel(Level.INFO);
		new MethylDbToFeatStatFiles().doMain(args);
	}

	public void doMain(String[] args)
	throws Exception
	{
		CmdLineParser parser = new CmdLineParser(this);
		// if you have a wider console, you could increase the value;
		// here 80 is also the default
		parser.setUsageWidth(80);
		try
		{
			parser.parseArgument(args);
		}
		catch (CmdLineException e)
		{
			System.err.println(e.getMessage());
			System.err.println(C_USAGE);
			// print the list of available options
			parser.printUsage(System.err);
			System.err.println();
			return;
		}

		List<String> tablePrefixes = new ArrayList<String>(3);
		List<String> featFns = new ArrayList<String>(100);
		boolean onFeats = false;
		for (String arg : arguments)
		{
			if (arg.equals(","))
			{
				//System.err.println("Found comma!");
				onFeats = true;
			}
			else
			{
				List<String> l = (onFeats) ? featFns : tablePrefixes;
				l.add(arg);
			}
		}
		
		if( (tablePrefixes.size()<1) || (featFns.size()<1))
		{
			System.err.println(C_USAGE);
			parser.printUsage(System.err);
			return;
		}
		
		int nS = tablePrefixes.size();
		int nFeatTypes = featFns.size();

	
		// Start the actual work.  Go through each GFF file.
		int onFeatType = 0;
		for (String featFn : featFns)
		{
			String featsFnBase = (new File(featFn)).getName();
			featsFnBase = featsFnBase.replaceFirst(".g[ft]f", "");
			Logger.getLogger(Logger.GLOBAL_LOGGER_NAME).info("On feature " + onFeatType + "/" + nFeatTypes + " " + featsFnBase);

			// Setup writer
			String outFn = String.format("%s.stats", featsFnBase);
			if (this.onechrom) outFn += ".onechrom";
			if (this.centeredSize>0) outFn += ".centered" + this.centeredSize;
			if (this.featShoreEnd>0) outFn += String.format(".featShoreStart%d-featShoreEnd%d", this.featShoreStart, this.featShoreEnd);
			ListUtils.setDelim("-");
			if (this.featFilters.size()>0) outFn += "." + ListUtils.excelLine(this.featFilters);
			outFn += ".csv";
			PrintWriter writer = new PrintWriter(new FileOutputStream(outFn));


			onFeatType++;
			ChromFeatures feats = new ChromFeatures(featFn, true);
			System.err.println("About to filter regions by size");
			if (this.maxFeatSize < Integer.MAX_VALUE) feats = feats.filterBySize(0, this.maxFeatSize);
			
			System.err.println("About to sort features");
			feats.sortFeatures();
			
			
			// The main work
			int featNum = 1;
			List<String> chroms = (this.onechrom) ? Arrays.asList("chr1") : MethylDbUtils.CHROMS;
			for (String chrStr : chroms)
			{
				int chrInt = feats.chrom_from_public_str(chrStr); 
				Iterator featit = feats.featureIterator(chrInt);
				System.err.println("Processing " + chrStr);
				
				// Keep the previous and next records to check for overlaps
				SimpleGFFRecord prevRec = (featit.hasNext()) ? (SimpleGFFRecord)featit.next() : null;
				SimpleGFFRecord thisRec = (featit.hasNext()) ? (SimpleGFFRecord)featit.next() : null;
				SimpleGFFRecord nextRec = (featit.hasNext()) ? (SimpleGFFRecord)featit.next() : null;
				// Process the first one 
				if (prevRec != null) processFeat(null, prevRec, thisRec, writer, tablePrefixes, chrStr, chrInt, feats);
				FEAT: do
				{
					//writer.println("On record " + featNum);
					try
					{
						if (thisRec != null) processFeat(prevRec, thisRec, nextRec, writer, tablePrefixes, chrStr, chrInt, feats);
						
						// Increment
						prevRec = thisRec;
						thisRec = nextRec;
						nextRec = (featit.hasNext()) ? (SimpleGFFRecord)featit.next() : null; 
					}
					catch (Exception e)
					{
						System.err.println("Couldn't complete feature " + featFn + ":\n" + e.toString());
						e.printStackTrace();
					}
					
					featNum++;

				} while (nextRec != null);

				// And the last one 
				if (thisRec != null) processFeat(prevRec, thisRec, null, writer, tablePrefixes, chrStr, chrInt, feats);
			
			}


			writer.close();

		}
	}

	protected void processFeat(GFFRecord prevRec, GFFRecord thisRec, GFFRecord nextRec, PrintWriter outPw, List<String> methylTables, String chrStr, int chrInt, ChromFeatures feats)
	throws Exception
	{
		
		List<String> cols = new ArrayList<String>(50);
		
		StrandedFeature.Strand featStrand = thisRec.getStrand();
		String featName = null; // rec.getSeqName();
		int nM = methylTables.size();


		int featS = thisRec.getStart();
		int featE = thisRec.getEnd();
		int featC = (int)Math.round((double)(featS+featE)/2.0);
		int flankS = featS;
		int flankE = featE;
		if (this.centeredSize > 0)
		{
			flankS = featC - this.centeredSize;
			flankE = featC + this.centeredSize;
		}
		
		if (this.featShoreEnd>0)
		{
			if (featStrand==StrandedFeature.NEGATIVE)
			{
				flankS = flankE + this.featShoreStart;
				flankE = flankE + this.featShoreEnd;
			}
			else
			{
				flankE = flankS - this.featShoreStart;
				flankS = flankS - this.featShoreEnd;
			}
			
			// Check if this overlaps another feature?
			SimpleGFFRecord shoreRec = new SimpleGFFRecord(thisRec);
			shoreRec.setStart(flankS);
			shoreRec.setEnd(flankE);
			if (((prevRec!=null) && MiscUtils.GffRecordsOverlap(shoreRec, prevRec)) ||
					(nextRec!=null) && MiscUtils.GffRecordsOverlap(shoreRec, nextRec))
			{
//				System.err.printf("Shore overlaps next feat, returning\n\tprev=%s\n\tthis=%s\n\tnext=%s\n",
//						GFFUtils.gffBetterString(prevRec),
//						GFFUtils.gffBetterString(thisRec),
//						GFFUtils.gffBetterString(nextRec));
				return;
			}
		}
		
		cols.add(Integer.toString(chrInt));
		cols.add(Integer.toString(flankS));
		cols.add(Integer.toString(flankE));

		// Expression
		double normalExp = Double.NaN, tumorExp = Double.NaN;
		normalExp = MethylDbUtils.fetchMeanExpression(chrStr, GFFUtils.getGffRecordName(thisRec), "N14838B");
		tumorExp = MethylDbUtils.fetchMeanExpression(chrStr, GFFUtils.getGffRecordName(thisRec), "(T14838A+T14838B)/2");
		cols.add(Double.toString(normalExp));
		cols.add(Double.toString(tumorExp));
		
		
		// Setup query params
		MethylDbQuerier params = new MethylDbQuerier();
		for (String featFilter : this.featFilters)
		{
			params.addFeatFilter(featFilter);
		}
		params.setMinCTreads(this.minCTreads);
		params.setMaxOppstrandAfrac(this.maxOppStrandAfrac);
		params.addRangeFilter(chrStr,flankS,flankE);

		// And stream through summarizers
		CpgIteratorMultisample cpgit = new CpgIteratorMultisample(params, methylTables);
		int numCpgs = cpgit.getCurNumRows();

		for (int m = 0; m < nM; m++)
		{
			// Initialize the first time
			if (this.methSummarizers[m] == null) 
			{
				this.methSummarizers[m] = new CpgMethLevelSummarizer(params);
				this.coverageSummarizers[m] = new CpgCoverageSummarizer(params);
			}
			this.methSummarizers[m].reset();
			this.coverageSummarizers[m].reset();
		}
		
		while (cpgit.hasNext()) 
		{
			Cpg[] cpgs = cpgit.next();

			for (int m = 0; m < nM; m++)
			{
				this.methSummarizers[m].streamCpg(cpgs[m]);
				this.coverageSummarizers[m].streamCpg(cpgs[m]);
			}
		}
		
		// Take averages
		for (int m = 0; m < nM; m++)
		{
			double meanM = this.methSummarizers[m].getValMean(true);
			double meanReads = this.coverageSummarizers[m].getValMean(true);
			cols.add(Integer.toString(numCpgs));
			cols.add(Double.toString(meanM));
			cols.add(Double.toString(meanReads));
		}
		
		
		
		ListUtils.setDelim(",");
		outPw.append( ListUtils.excelLine(cols));
		outPw.println();
	}
	

	
	
}
