package edu.usc.epigenome.scripts;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;

import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;

import edu.usc.epigenome.genomeLibs.MethylDb.Cpg;
import edu.usc.epigenome.genomeLibs.MethylDb.CpgIterator;
import edu.usc.epigenome.genomeLibs.MethylDb.MethylDbQuerier;



public class MethylDbToWig_DEPR {

	private static final String C_USAGE = "Use: MethylDbToWig -withinFeat featType -tablePrefix " + MethylDbQuerier.DEFAULT_METHYL_TABLE_PREFIX + 
	" -minCTreads 10 -maxOppStrandAfrac 0.10 -noNonconvFilter chr [startPos] [endPos]";
	
    @Option(name="-noNonconvFilter",usage="override the nonconversion filter (default false)")
    protected boolean noNonconvFilter = false;
    @Option(name="-tablePrefix",usage="Prefix for DB table (default " + MethylDbQuerier.DEFAULT_METHYL_TABLE_PREFIX + ")")
    protected String tablePrefix = null;
    @Option(name="-withinFeat",usage="A featType from the features table")
    protected String withinFeat = null;
    @Option(name="-minCTreads",usage="Minimum number of C or T reads to count as a methylation value")
    protected int minCTreads = 0;
    @Option(name="-maxOppStrandAfrac",usage="As on the opposite strand are evidence for mutation or SNP. " +
    		"This sets a maximum number of observed As on the opposite strand (default 0.1)")
    protected double maxOppStrandAfrac = 0.1;
	// receives other command line parameters than options
	@Argument
	private List<String> arguments = new ArrayList<String>();
	
	public static void main(String[] args)
	throws Exception	
	{
		new MethylDbToWig_DEPR().doMain(args);
	}

	public void doMain(String[] args)
	throws Exception
	{
		CmdLineParser parser = new CmdLineParser(this);
		// if you have a wider console, you could increase the value;
		// here 80 is also the default
		int chrSt = 0, chrEnd = 0;
		String chr;
		parser.setUsageWidth(80);
		try
		{
			parser.parseArgument(args);

			if(arguments.size() < 1 ) {
				System.err.println(C_USAGE);
				System.exit(1);
			}

			chr = arguments.get(0);
			if (!chr.startsWith("chr")) chr = "chr" + chr;
			if (arguments.size()>1)
			{
				chrSt = Integer.parseInt(arguments.get(1));
				chrEnd = Integer.parseInt(arguments.get(2));
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
		

		Logger.getLogger(Logger.GLOBAL_LOGGER_NAME).setLevel(Level.SEVERE);
		
		
		MethylDbQuerier params = new MethylDbQuerier();
		if (this.tablePrefix != null) params.methylTablePrefix = this.tablePrefix;
		params.setMinCTreads(this.minCTreads);
		params.setUseNonconversionFilter(!this.noNonconvFilter);
		params.setMaxOppstrandAfrac(this.maxOppStrandAfrac);
		if (this.withinFeat!=null) params.addFeatFilter(this.withinFeat, 0);
		
		if (arguments.size()>1)
		{
			params.addRangeFilter(chr, chrSt, chrEnd);
		}
		else
		{
			params.addRangeFilter(chr);
		}

		Iterator<Cpg> it;
		it = new CpgIterator(params);
		
		int count = 0;
		System.out.printf("track type=wiggle_0 name=%s description=%s\n", "test", "test");
		System.out.printf("variableStep chrom=%s\n",chr);
		
		while (it.hasNext())
		{
			Cpg cpg = it.next();
			if (!Double.isNaN(cpg.fracMeth(!this.noNonconvFilter)))
			{
				System.out.println(cpg.variableStepWigLine(!this.noNonconvFilter));
				//System.err.println(cpg.toString());
			}
			count++;
		}
		
		System.err.printf("Found %d items\n",count);
	
	}

}
