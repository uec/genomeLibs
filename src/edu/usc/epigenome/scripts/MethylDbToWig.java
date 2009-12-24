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



public class MethylDbToWig {

	private static final String C_USAGE = "Use: MethylDbToWig -noNonconvFilter chr [startPos] [endPos]";
	
    @Option(name="-noNonconvFilter",usage="override the nonconversion filter (default false)")
    protected boolean noNonconvFilter = false;
	// receives other command line parameters than options
	@Argument
	private List<String> arguments = new ArrayList<String>();
	
	public static void main(String[] args)
	throws Exception	
	{
		new MethylDbToWig().doMain(args);
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

			chr = arguments.get(0);
			if (arguments.size()>1)
			{
				chrSt = Integer.parseInt(arguments.get(1));
				chrEnd = Integer.parseInt(arguments.get(2));
			}

			if(arguments.size() < 1 ) {
				System.err.println(C_USAGE);
				System.exit(1);
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
		
		
		Iterator<Cpg> it;
		
		if (arguments.size()>1)
		{
			it = new CpgIterator(chr, chrSt, chrEnd);
		}
		else
		{
			it = new CpgIterator(chr, chrSt, chrEnd);
		}
		
		int count = 0;
		System.out.printf("track type=wiggle_0 name=%s description=%s\n", "test", "test");
		System.out.printf("variableStep chrom=%s\n",chr);
		
		while (it.hasNext())
		{
			Cpg cpg = it.next();
			if (!Double.isNaN(cpg.fracMeth(!this.noNonconvFilter)))
			{
				System.out.println(cpg.variableStepWigLine(!this.noNonconvFilter));
			}
			count++;
		}
		
		System.err.printf("Found %d items\n",count);
	
	}

}
