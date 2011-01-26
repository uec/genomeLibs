package edu.usc.epigenome.scripts;

import java.io.FileOutputStream;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;

import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;

import edu.usc.epigenome.genomeLibs.GoldAssembly;
import edu.usc.epigenome.genomeLibs.MatUtils;
import edu.usc.epigenome.genomeLibs.MethylDb.Cpg;
import edu.usc.epigenome.genomeLibs.MethylDb.CpgIterator;
import edu.usc.epigenome.genomeLibs.MethylDb.CpgIteratorMultisample;
import edu.usc.epigenome.genomeLibs.MethylDb.MethylDbQuerier;
import edu.usc.epigenome.genomeLibs.MethylDb.MethylDbUtils;



public class MethylDbToBareWigs {

	private static final String C_USAGE = "Use: MethylDbToWigs -outPrefix -withinFeat featType -table methylCGsRich_normal010310_ " + 
	" -table methylCGsRich_tumor011010_ -minCTreads 2 -maxOppStrandAfrac 0.10 -noNonconvFilter chr [startPos] [endPos]";

	@Option(name="-noNonconvFilter",usage="override the nonconversion filter (default false)")
	protected boolean noNonconvFilter = false;
	@Option(name="-withinFeat",usage="A featType from the features table")
	protected String withinFeat = null;
	@Option(name="-table",usage="A table prefix")
	protected List<String> tables = new ArrayList<String>(10);
	@Option(name="-outPrefix",usage="Output files will have this name")
	protected String outPrefix = "wiggleTester";
	@Option(name="-minCTreads",usage="Minimum number of C or T reads to count as a methylation value")
	protected int minCTreads = 2;
	@Option(name="-extraCols",usage="print C and T read counts")
	protected boolean extraCols = false;
	@Option(name="-maxOppStrandAfrac",usage="As on the opposite strand are evidence for mutation or SNP. " +
	"This sets a maximum number of observed As on the opposite strand (default 0.1)")
	protected double maxOppStrandAfrac = 0.1;
	// receives other command line parameters than options
	@Argument
	private List<String> arguments = new ArrayList<String>();

	public static void main(String[] args)
	throws Exception	
	{
		new MethylDbToBareWigs().doMain(args);
	}

	public void doMain(String[] args)
	throws Exception
	{
		boolean chrEndSet = false;
		CmdLineParser parser = new CmdLineParser(this);
		// if you have a wider console, you could increase the value;
		// here 80 is also the default
		int chrSt = 0, chrEnd = GoldAssembly.chromLengthStatic("chr1", "hg18");
		List<String> chrs;
		parser.setUsageWidth(80);
		try
		{
			parser.parseArgument(args);

			if(arguments.size() < 1 ) {
				chrs = MethylDbUtils.CHROMS;
			}
			else
			{
				String chr = arguments.get(0);
				if (!chr.startsWith("chr")) chr = "chr" + chr;
				chrs = new ArrayList<String>(1);
				chrs.add(chr);
				if (arguments.size()>1)
				{
					chrSt = Integer.parseInt(arguments.get(1));
					chrEnd = Integer.parseInt(arguments.get(2));
					chrEndSet = true;
				}
			}

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

		if (chrSt <= 0) chrSt = 1;
		
		if (tables.size()==0){
			System.err.println("Error, must supply at least one table");
			System.err.println(C_USAGE);
			// print the list of available options
			parser.printUsage(System.err);
			System.err.println();
			return;
		}

		Logger.getLogger(Logger.GLOBAL_LOGGER_NAME).setLevel(Level.SEVERE);


		MethylDbQuerier params = new MethylDbQuerier();
		params.setMinCTreads(this.minCTreads);
		params.setUseNonconversionFilter(!this.noNonconvFilter);
		params.setMaxOppstrandAfrac(this.maxOppStrandAfrac);
		if (this.withinFeat!=null) params.addFeatFilter(this.withinFeat, 0);

		if (chrs.size()==1) this.outPrefix += "-" + chrs.get(0);

		int nTabs = tables.size();
		PrintWriter[] pws = new PrintWriter[nTabs];
		for (int i = 0; i < nTabs; i++)
		{
			String tab = tables.get(i);

			// Setup output files.
			String outfn = String.format("%s.%s.minCTreads%d.maxOppAfrac%.3f.bare.wig", 
					outPrefix,tab,minCTreads,this.maxOppStrandAfrac);
			PrintWriter pw = new PrintWriter(new FileOutputStream(outfn));
			pws[i] = pw;
			// Use -2 for view limit so you can see ones at 0
			pw.printf("track type=wiggle_0 name=%sbare description=%sbare color=204,102,0 visibility=full " +
					" graphType=points autoScale=off alwaysZero=off maxHeightPixels=64:32:10 viewLimits=-2:100\n", tab, tab);
		}


		for (String chr : chrs)
		{

			if (!chrEndSet)
			{
				chrEnd = GoldAssembly.chromLengthStatic(chr, "hg18");
			}

			for (int i = 0; i < nTabs; i++)
			{
				String rawspanline = String.format("variableStep chrom=%s\n",chr);
				pws[i].append(rawspanline);
			}

			// If memory gets to be a problem with whole genome loads, try this trick to walk the chromosome.		
			// Stupid JDBC tries to load entire chromosome into memory at once,
			// which is too much.
			final int STEP = (int)1E7;
			for (int c = chrSt; c < (chrEnd+STEP); c += STEP)
			{
				Logger.getLogger(Logger.GLOBAL_LOGGER_NAME).severe("Starting chrom " + chr + " from " + c+1 + " to " + c+STEP + "\n");

				params.clearRangeFilters();
				params.addRangeFilter(chr, c+1, Math.min(chrEnd, c+STEP));					
				CpgIteratorMultisample cpgit = new CpgIteratorMultisample(params, this.tables);
				int numCpgs = cpgit.getCurNumRows();

				CPG: for (int i = 0; i < numCpgs; i++)
				{
					Cpg[] cpgs = null;
					try
					{
						cpgs = cpgit.next();
					}
					catch (Exception e)
					{
						System.err.printf("cpgit.getCurNumRows() said we had %d rows, but i am failing on the %dth!\n%s",
								numCpgs, i, e.toString());
						continue CPG;
					}

					for (int j = 0; j < nTabs; j++)
					{
						double meth = cpgs[j].fracMeth(!this.noNonconvFilter);
						if (!Double.isNaN(meth))
						{					
							int pos = cpgs[j].chromPos;
							pws[j].printf("%d\t%d",pos,(int)Math.round(100.0*meth));
							if (extraCols) pws[j].printf("\t%d\t%d",cpgs[j].cReads,cpgs[j].cReads+cpgs[j].tReads);
							pws[j].println();
						}
					}

				} // Cpg
			} // Segment
		} // Chr

		for (int i = 0; i < nTabs; i++)
		{
			pws[i].close();
		}		

	} // Main

	public static double safeVal(double v)
	{
		return ((Double.isNaN(v) || Double.isInfinite(v)) ? 0.0 : v);
	}

}





