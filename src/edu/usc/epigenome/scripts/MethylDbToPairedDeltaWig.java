package edu.usc.epigenome.scripts;

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

import com.sun.org.apache.xalan.internal.xsltc.compiler.Pattern;

import edu.usc.epigenome.genomeLibs.MethylDb.Cpg;
import edu.usc.epigenome.genomeLibs.MethylDb.CpgIterator;
import edu.usc.epigenome.genomeLibs.MethylDb.CpgIteratorMultisample;
import edu.usc.epigenome.genomeLibs.MethylDb.MethylDbQuerier;



public class MethylDbToPairedDeltaWig {

	private static final String C_USAGE = "Use: MethylDbToPairedDeltaWig " +  
	" -minCTreads 10 -maxOppStrandAfrac 0.10 -noNonconvFilter chr [startPos] [endPos] table1prefix table2prefix";
	
    @Option(name="-noNonconvFilter",usage="override the nonconversion filter (default false)")
    protected boolean noNonconvFilter = false;
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
		new MethylDbToPairedDeltaWig().doMain(args);
	}

	public void doMain(String[] args)
	throws Exception
	{
		CmdLineParser parser = new CmdLineParser(this);
		// if you have a wider console, you could increase the value;
		// here 80 is also the default
		int chrSt = -1, chrEnd = -1;
		String chr, table1prefix, table2prefix;
		parser.setUsageWidth(80);
		try
		{
			parser.parseArgument(args);

			if(arguments.size() < 1 ) {
				System.err.println(C_USAGE);
				System.exit(1);
			}

			chr = arguments.remove(0);
			if (!chr.startsWith("chr")) chr = "chr" + chr;
			if (arguments.size()==4)
			{
				chrSt = Integer.parseInt(arguments.remove(0));
				chrEnd = Integer.parseInt(arguments.remove(0));
			}
			table1prefix = arguments.remove(0);
			table2prefix = arguments.remove(0);
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
		params.minCTreads = this.minCTreads;
		params.useNonconversionFilter = !this.noNonconvFilter;
		params.maxOppstrandAfrac = this.maxOppStrandAfrac;
		if (chrSt >= 0)
		{
			params.addRangeFilter(chr, chrSt, chrEnd);
		}
		else
		{
			params.addRangeFilter(chr);
		}

		Iterator<Cpg[]> it;
		it = new CpgIteratorMultisample(params, Arrays.asList(table1prefix,table2prefix));
		
		int count = 0;
		System.out.printf("track type=wiggle_0 name=%s description=%s\n", "test", "test");
		System.out.printf("variableStep chrom=%s\n",chr);
		
		while (it.hasNext())
		{
			Cpg[] cpgs = it.next();
			int ns = cpgs.length;
			
			
			double a = cpgs[0].fracMeth(!this.noNonconvFilter);
			double b = cpgs[1].fracMeth(!this.noNonconvFilter);
			double delta = b - a;
			if (!Double.isNaN(a) && !Double.isNaN(b))
			{
				System.out.println(cpgs[0].chromPos + "\t" + delta);
			}
			
//			for (int i = 0; i < ns; i++)
//			{
//				Cpg cpg = cpgs[i];
//				if (i>0) System.out.print("\t");
//			
//				System.out.print(cpg.variableStepWigLine(!this.noNonconvFilter));
//				// System.out.print(cpg.toString());
//			}
//			System.out.println();

			count++;
		}
		
		System.err.printf("Found %d items\n",count);
	
	}

}
