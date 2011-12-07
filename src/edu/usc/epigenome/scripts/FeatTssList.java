package edu.usc.epigenome.scripts;

import java.io.File;
import java.io.FileOutputStream;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;
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
import edu.usc.epigenome.genomeLibs.MethylDb.Cpg;
import edu.usc.epigenome.genomeLibs.MethylDb.CpgIterator;
import edu.usc.epigenome.genomeLibs.MethylDb.CpgIteratorMultisample;
import edu.usc.epigenome.genomeLibs.MethylDb.MethylDbQuerier;
import edu.usc.epigenome.genomeLibs.MethylDb.MethylDbUtils;


public class FeatTssList {

	private static final String C_USAGE = "Use: FeatTssList";
	
//	@Option(name="-skipUnoriented",usage="If set, skip any unoriented feature (default false)")
//	protected boolean skipUnoriented = false;
	@Option(name="-chrom",multiValued=true,usage="One or more chroms, eg. -chrom chr1 -chrom chr5")
	protected List<String> chrs = new ArrayList<String>(25);
	@Option(name="-feature",multiValued=true,usage="One or more features from features_ tables")
	protected List<String> features = new ArrayList<String>(25);
	@Option(name="-flank",multiValued=false,usage="Flanking sequence to search for each feature (default 0=exact overlap)")
	protected int flank = 0;


	// receives other command line parameters than options
	@Argument
	private List<String> arguments = new ArrayList<String>();


	// class vars


	/**
	 * @param args
	 */
	public static void main(String[] args)
	throws Exception	
	{
		Logger.getLogger(Logger.GLOBAL_LOGGER_NAME).setLevel(Level.INFO);
		new FeatTssList().doMain(args);
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

		if (arguments.size()>0)
		{
			System.err.println(C_USAGE);
			parser.printUsage(System.err);
			return;
		}
		
		if (features.size()<1)
		{
			System.err.println("Must supply at least 1 -feature argument");
			System.err.println(C_USAGE);
			parser.printUsage(System.err);
			return;
		}		
		
		if (chrs.size()==0) chrs = MethylDbUtils.CHROMS;
		


		// Go through target feats one by one

		for (String chrStr : chrs)
		{
			// Some of my files have these forms.
			if (chrStr.equalsIgnoreCase("chrX")) chrStr = "chrX";
			if (chrStr.equalsIgnoreCase("chry")) chrStr = "chrY";
			if (chrStr.equalsIgnoreCase("chrm")) chrStr = "chrM";

			
			FeatDbQuerier tssParams = new FeatDbQuerier();
			//tssParams.outputSymbolInsteadOfRefseq = true;
			tssParams.addRangeFilter(chrStr, 0, 0); // Set to 0,0 for whole chromosome
			tssParams.addFeatFilter("tss");
			FeatIterator tssFeats = new FeatIterator(tssParams);
			while (tssFeats.hasNext())
			{
				GFFRecord tssRec = tssFeats.next();
				String tssName = GFFUtils.getGffRecordName(tssRec);		
				

				FeatDbQuerier params = new FeatDbQuerier();
				params.addRangeFilter(chrStr, tssRec.getStart()-this.flank, tssRec.getEnd()+this.flank);
				for (int n = 0; n < features.size(); n++)
				{
					String featType = features.get(n);
					params.addFeatFilter(featType);
				}
				FeatIterator feats = new FeatIterator(params);
				if (feats.hasNext())
				{
					System.out.println(tssName);
				}

				
			}
		}
	}
	
}
