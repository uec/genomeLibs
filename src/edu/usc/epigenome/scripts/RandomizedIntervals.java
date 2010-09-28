package edu.usc.epigenome.scripts;

import java.io.File;
import java.io.FileOutputStream;
import java.io.PrintWriter;
import java.util.ArrayList;
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

import sun.tools.tree.ThisExpression;

import edu.usc.epigenome.genomeLibs.GFFUtils;
import edu.usc.epigenome.genomeLibs.GoldAssembly;
import edu.usc.epigenome.genomeLibs.ListUtils;
import edu.usc.epigenome.genomeLibs.MatUtils;
import edu.usc.epigenome.genomeLibs.ChromScores.ChromScoresArrayInt;
import edu.usc.epigenome.genomeLibs.FeatAligners.FeatAligner;
import edu.usc.epigenome.genomeLibs.FeatAligners.FeatAlignerAveraging;
import edu.usc.epigenome.genomeLibs.FeatAligners.FeatAlignerEachfeat;
import edu.usc.epigenome.genomeLibs.FeatDb.FeatDbQuerier;
import edu.usc.epigenome.genomeLibs.FeatDb.FeatIterator;
import edu.usc.epigenome.genomeLibs.GenomicRange.GenomicRange;
import edu.usc.epigenome.genomeLibs.GenomicRange.GenomicRangeRandomizer;
import edu.usc.epigenome.genomeLibs.MethylDb.Cpg;
import edu.usc.epigenome.genomeLibs.MethylDb.CpgIterator;
import edu.usc.epigenome.genomeLibs.MethylDb.CpgIteratorMultisample;
import edu.usc.epigenome.genomeLibs.MethylDb.MethylDbQuerier;
import edu.usc.epigenome.genomeLibs.MethylDb.MethylDbUtils;


public class RandomizedIntervals {

	private static final String C_USAGE = "Use: RandomizedIntervals -numTrials 1 -matchDistanceGtf distmatch.gtf -maskInGtf maskIn.gtf -maskOutGtf maskOutGft -genome hg18 mimic1.gtf mimic2.gtf ...";
	
//	@Option(name="-skipUnoriented",usage="If set, skip any unoriented feature (default false)")
//	protected boolean skipUnoriented = false;
	@Option(name="-mimicGtf",multiValued=false,usage="Flanking sequence for each feature")
	protected String mimicGtf = null;
	@Option(name="-maskInGtf",multiValued=false,usage="Include only features within these regions (default off)")
	protected String maskInGtf = null;
	@Option(name="-maskOutGtf",multiValued=false,usage="Include only features within these regions (default off)")
	protected String maskOutGtf = null;
	@Option(name="-matchDistanceGtf",multiValued=false,usage="Match distances, for instance distances to TSS is useful (default off)")
	protected String matchDistanceGtf = null;
	@Option(name="-genome",multiValued=false,usage="Genome assembly (default hg18)")
	protected String genome = "hg18";
	@Option(name="-numTrials",multiValued=false,usage="If greater than 1, we do multiple independent trials for each mimic file")
	protected int numTrials = 1;

	// receives other command line parameters than options
	@Argument
	private List<String> mimicFns = new ArrayList<String>();


	// class vars


	/**
	 * @param args
	 */
	public static void main(String[] args)
	throws Exception	
	{
		Logger.getLogger(Logger.GLOBAL_LOGGER_NAME).setLevel(Level.INFO);
		new RandomizedIntervals().doMain(args);
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

		if (mimicFns.size()<1)
		{
			System.err.println("Must provide at least one mimic file");
			System.err.println(C_USAGE);
			parser.printUsage(System.err);
			return;
		}
		
		
		GenomicRangeRandomizer randomizer = new GenomicRangeRandomizer(this.genome, this.maskInGtf, this.maskOutGtf, this.matchDistanceGtf);

		for (String mimicFn : mimicFns)
		{
			ChromFeatures cf = new ChromFeatures(mimicFn, true);
			
			
			for (int t = 1; t <= this.numTrials; t++)
			{
				String outFn = mimicFn;
				outFn = outFn.replaceAll(".g[tf]f$", "");
				outFn += ".randomizedLocs.trial" + t + ".gtf";
				PrintWriter outwriter = new PrintWriter(new FileOutputStream(outFn));

				
				System.err.printf("On trial %d\n", t);
				Iterator featIt = cf.featureIterator();
				while (featIt.hasNext())
				{
					SimpleGFFRecord target = (SimpleGFFRecord) featIt.next();
					int len = target.getEnd()-target.getStart()+1;
					GenomicRange randRange = randomizer.generateRandomRange(len);
					randRange.setScore(target.getScore());
					randRange.setStrand(target.getStrand());
					outwriter.println(randRange.gffStr());
				}
				
				outwriter.close();

			}
		}
		
	}
	
	
}
