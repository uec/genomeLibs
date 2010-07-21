package edu.usc.epigenome.scripts;

import java.io.File;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;

import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;

import com.sun.org.apache.xalan.internal.xsltc.compiler.Pattern;

import edu.usc.epigenome.genomeLibs.GoldAssembly;

import edu.usc.epigenome.genomeLibs.MethylDb.Cytosine;
import edu.usc.epigenome.genomeLibs.MethylDb.CytosineIterator;
import edu.usc.epigenome.genomeLibs.MethylDb.MethylDbQuerier;

public class MethylDbToGpcDensities {

	private static final String C_USAGE = "Use: MethylDbToGpcDensities -tablePrefix " + MethylDbQuerier.DEFAULT_METHYL_TABLE_PREFIX + 
	" stepSize windSize chr [startPos] [endPos]";
	
    @Option(name="-tablePrefix",usage="Prefix for DB table (default " + MethylDbQuerier.DEFAULT_METHYL_TABLE_PREFIX + ")")
    //protected String tablePrefix = null;
    protected String tablePrefix = "GpC_test_";
	// receives other command line parameters than options
	@Argument
	private List<String> arguments = new ArrayList<String>();
	
	public static void main(String[] args) 
	throws Exception{
		// TODO Auto-generated method stub
		new MethylDbToGpcDensities().doMain(args);
	}

	public void doMain(String[] args)
	throws Exception
	{
		CmdLineParser parser = new CmdLineParser(this);
		// if you have a wider console, you could increase the value;
		// here 80 is also the default
		int windSize = 0, stepSize = 0, chrStart = -1, chrEnd = -1;
		String chr;
		parser.setUsageWidth(80);
		try
		{
			parser.parseArgument(args);

			if(arguments.size() < 3 ) {
				System.err.println(C_USAGE);
				System.exit(1);
			}

			stepSize = Integer.parseInt(arguments.get(0));
			windSize = Integer.parseInt(arguments.get(1));
			chr = arguments.get(2);
			if (!chr.startsWith("chr")) chr = "chr" + chr;
			if (arguments.size()>3)
			{
				chrStart = Integer.parseInt(arguments.get(3));
				chrEnd = Integer.parseInt(arguments.get(4));
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
		

		Logger.getLogger(Logger.GLOBAL_LOGGER_NAME).setLevel(Level.INFO);
		
		MethylDbQuerier params = new MethylDbQuerier();
		if (this.tablePrefix != null) params.methylTablePrefix = this.tablePrefix;
		params.setMinCTreads(0);
		params.setUseNonconversionFilter(false);

		if (chrStart<0)
		{
			chrStart = 1;
			chrEnd = GoldAssembly.chromLengthStatic(chr, "hg18");
			//System.out.println(chrEnd);
		}
		
		if (windSize<=0) windSize = 1000;
		if (stepSize<=0) stepSize = 1000;

		int halfWind = (int)Math.floor((double)windSize/(double)2);
		
		System.out.printf("track type=wiggle_0 name=%s description=%s\n", "GpcDensity", "GpcDensity");
		System.out.printf("fixedStep chrom=%s start=%d step=%d span=%d\n",chr,chrStart,stepSize,windSize);

		int numQueries = 0;
		String fn = tablePrefix + chr +"_density.txt";
		PrintWriter writer = new PrintWriter(new File(fn));
		for (int c = chrStart; c < chrEnd; c+=stepSize)
		{
			if ((++numQueries % 10000)==0) Logger.getLogger(Logger.GLOBAL_LOGGER_NAME).info("On query " + numQueries);
			
			
			int windS = c - halfWind;
			int windE = c + halfWind - 1;
			if (windS<chrStart) windS = chrStart;
			if (windE>chrEnd) windE = chrEnd;
			int windLen = windE-windS+1;
			Logger.getLogger(Logger.GLOBAL_LOGGER_NAME).fine(String.format("%d,%d,%d\t%d\n", windS,c,windE,windLen));
			System.out.printf("%s\t%d\t%d\n",chr, windS, windE);
			params.clearRangeFilters();
			params.addRangeFilter(chr, windS, windE);
			
			CytosineIterator it = new CytosineIterator(params);
			int numGpc = 0;
			
			//while ((it.getString("preBaseRefUpperCase") == "G") && (it.getString("nextBaseRefUpperCase") != "G")) {
				
			
			//int numGpcs = it.getCurNumRows();
			//double numGpcs = (double)numCs/2.0;  // We have 2 in the database for each Gpc (one per strand)
			
			//float methyDens =  it.getMethyDens();
			double methySum = 0;
			while(it.hasNext()){
				Cytosine methyValue = it.next();
				if((methyValue.getPreBaseRef() == 'G') && (methyValue.getNextBaseRef() != 'G')){
					
					methySum += methyValue.fracMeth();
					numGpc++;
				}
			}
			int gpcDens = (int)Math.round(100.0 * numGpc / ((double)windLen *2));  // Make it percent of sequence
			double methyDens = methySum / (numGpc*100);
			//System.out.printf("%d\t%d\t%d\t%f\n",gpcDens, numGpcs, windLen, methyDens);
//			System.out.println(dens);
			//String fn = prefix + sampleName + "_" + chr + ".txt";
			
			
			writer.printf("%d\t%f\t%d\t%d\n", windS, methyDens, numGpc, gpcDens);
			
		}
		writer.close();

	}
}
