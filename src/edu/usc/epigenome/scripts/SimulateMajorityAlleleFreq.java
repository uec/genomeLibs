package edu.usc.epigenome.scripts;

import java.util.ArrayList;
import java.util.List;

import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;

import edu.usc.epigenome.genomeLibs.MatUtils;



public class SimulateMajorityAlleleFreq {

	final private static String USAGE = "SimulateMajorityAlleleFreq [opts] > means.out";

	/**
	 * object vars
	 */

	
	/**
	 * @param args
	 */

	@Option(name="-numIterations",usage="run the simulation this many times (default 1)")
	protected int numIterations = 1;
	@Option(name="-numSamples",usage="number of samples to average for each iteration (default 10000)")
	protected int numSamples = 10000;
	@Option(name="-allelesPerSample",usage="number of alleles per sample (default 4)")
	protected int allelesPerSample = 4;
	@Option(name="-alleleFreq",usage="Frequency of allele A (default 0.5)")
	protected double alleleFreq = 0.5;
	@Option(name="-debug",usage=" Debugging statements (default false)")
	protected boolean debug = false;

	
	// receives other command line parameters than options
	@Argument
	private List<String> stringArgs = new ArrayList<String>();

	
	
	/**
	 * @param args
	 */
	public static void main(String[] args)
	throws Exception
	{
		new SimulateMajorityAlleleFreq().doMain(args);
	}

	public void doMain(String[] args)
	throws Exception {

		CmdLineParser parser = new CmdLineParser(this);
		// if you have a wider console, you could increase the value;
		// here 80 is also the default
		parser.setUsageWidth(80);
		try
		{
			parser.parseArgument(args);
			if (stringArgs.size() > 0) throw new CmdLineException(USAGE);
			
		}
		catch (CmdLineException e)
		{
			System.err.println(e.getMessage());
			// print the list of available options
			parser.printUsage(System.err);
			System.err.println();
			return;
		}
		
		
		double[] means = new double[this.numIterations]; 

		for (int it_i = 0; it_i < this.numIterations; it_i++)
		{
			double fracTotal = 0.0;
			
			for (int i = 0; i < this.numSamples; i++)
			{
				int alleleAcount = 0;
				for (int a = 0; a < this.allelesPerSample; a++)
				{
					if (Math.random() > this.alleleFreq) alleleAcount++;
				}
				int alleleBcount = this.allelesPerSample-alleleAcount;
				int majorityAlleleCount = (alleleAcount>alleleBcount) ? alleleAcount : alleleBcount;
				
				double frac = (double)majorityAlleleCount/(double)this.allelesPerSample;
				fracTotal += frac;
			}
			
			double mean = (fracTotal / (double)this.numSamples);
			//System.err.format("%.3f\n", mean);
			means[it_i] = mean;
		}
		
		double meanOfMeans = MatUtils.nanMean(means);
		double stdevOfMeans = MatUtils.nanStdev(means);
		System.out.format("%.5f\t%.5f\n", meanOfMeans,stdevOfMeans);
		System.err.format("%.5f\t%.5f\n", meanOfMeans,stdevOfMeans);

		
	}


	
}
