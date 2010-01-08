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
import edu.usc.epigenome.genomeLibs.MatUtils;
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


public class FeatCounter {

	private static final String C_USAGE = "Use: FeatCounter larger.gtf smaller";
	
//	@Option(name="-skipUnoriented",usage="If set, skip any unoriented feature (default false)")
//	protected boolean skipUnoriented = false;
	@Option(name="-chrom",multiValued=true,usage="One or more chroms, eg. -chrom chr1 -chrom chr5")
	protected List<String> chrs = new ArrayList<String>(25);

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
		new FeatCounter().doMain(args);
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

		if (arguments.size()!=2)
		{
			System.err.println(C_USAGE);
			parser.printUsage(System.err);
			return;
		}
		
		if (chrs.size()==0) chrs = MethylDbUtils.CHROMS;
		String targetFn = arguments.get(0);
		String featType = arguments.get(1);


		ChromFeatures targetFeats = new ChromFeatures(targetFn, true);
		// Go through target feats one by one
		for (String chrStr : chrs)
		{
			// Some of my files have these forms.
			if (chrStr.equalsIgnoreCase("chrX")) chrStr = "chrX";
			if (chrStr.equalsIgnoreCase("chry")) chrStr = "chrY";
			if (chrStr.equalsIgnoreCase("chrm")) chrStr = "chrM";

			int nF = 0;
			Iterator targetit = targetFeats.featureIterator((new ChromFeatures()).chrom_from_public_str(chrStr));
			while (targetit.hasNext())
			{

				SimpleGFFRecord target = (SimpleGFFRecord) targetit.next();



				// Get the overlapping feats
				FeatDbQuerier params = new FeatDbQuerier();
				params.addRangeFilter(target.getSeqName(), target.getStart(), target.getEnd());
				params.addFeatFilter(featType);
				FeatIterator feats = new FeatIterator(params);

				while (feats.hasNext())
				{
					feats.next();
					nF++;
				}
			}
			System.out.printf("%s,%s,%d\n", featType, chrStr, nF);
		}

	}
	
	
}
