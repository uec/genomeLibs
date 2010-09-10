package edu.usc.epigenome.scripts;

import java.io.FileOutputStream;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;

import org.biojava.bio.seq.StrandedFeature;
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



public class MethylDbToWigs {

	private static final String C_USAGE = "Use: MethylDbToWigs -outPrefix -withinFeat featType -table1 methylCGsRich_normal010310_ " + 
		" -stepSize 50 -windSize 100 -minCpgs 6 -maxWindStretch 5000 -table2 methylCGsRich_tumor011010_ -minCTreads 2 -maxOppStrandAfrac 0.10 -noNonconvFilter chr [startPos] [endPos]";
	
    @Option(name="-bare",usage="Create the large bare wiggle files")
    protected boolean bare = false;
    @Option(name="-noNonconvFilter",usage="override the nonconversion filter (default false)")
    protected boolean noNonconvFilter = false;
    @Option(name="-table1",usage="Prefix for DB table (default methylCGsRich_normalM030510_)")
    protected String table1 = "methylCGsRich_normalM030510_";
    @Option(name="-table2",usage="Prefix for DB table (default methylCGsRich_tumorM030510_)")
    protected String table2 = "methylCGsRich_tumorM030510_";
    @Option(name="-withinFeat",usage="A featType from the features table")
    protected String withinFeat = null;
    @Option(name="-outPrefix",usage="Output files will have this name")
    protected String outPrefix = null;
    @Option(name="-stepSize",usage="step size for wiggle (50)")
    protected int stepSize = 50;
    @Option(name="-windSize",usage="starting window size (500)")
    protected int windSize = 500;
    @Option(name="-minCpgs",usage="minimum Cpgs in window (20)")
    protected int minCpgs = 20;
    @Option(name="-maxWindStretch",usage="maximum amount to stretch window to find minCpgs cpgs (5000)")
    protected int maxWindStretch = 10000;
    @Option(name="-minCTreads",usage="Minimum number of C or T reads to count as a methylation value")
    protected int minCTreads = 2;
    @Option(name="-maxOppStrandAfrac",usage="As on the opposite strand are evidence for mutation or SNP. " +
    		"This sets a maximum number of observed As on the opposite strand (default 0.1)")
    protected double maxOppStrandAfrac = 0.1;
	// receives other command line parameters than options
	@Argument
	private List<String> arguments = new ArrayList<String>();
	
	public static void main(String[] args)
	throws Exception	
	{
		new MethylDbToWigs().doMain(args);
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

			chrs = MethylDbUtils.TEST_CHROMS; // MethylDbUtils.CHROMS;
			if(this.outPrefix==null)
			{
				System.err.println(C_USAGE);
				// print the list of available options
				parser.printUsage(System.err);
				System.exit(0);
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
		
		Logger.getLogger(Logger.GLOBAL_LOGGER_NAME).setLevel(Level.SEVERE);
		
		
		MethylDbQuerier params = new MethylDbQuerier();
		params.setMinCTreads(this.minCTreads);
		params.setUseNonconversionFilter(!this.noNonconvFilter);
		params.setMaxOppstrandAfrac(this.maxOppStrandAfrac);
		if (this.withinFeat!=null) params.addFeatFilter(this.withinFeat, 0);
		
		if (chrs.size()==1) this.outPrefix += "-" + chrs.get(0);

		// Setup output files.
		List<PrintWriter> openFiles = new ArrayList<PrintWriter>(100);
		
		// Special cases for coloring
		String colorA = "204,153,102";
		String colorB = "204,102,0";
		if (this.table1.toUpperCase().contains("H1"))
		{
			colorA = "0,102,0";
			colorB = "0,102,255";
		}
		else if (this.table1.toUpperCase().contains("IMR90"))
		{
			colorB = "0,102,0";
			colorA = "0,102,255";			
		}

		String outfnAs = this.outPrefix + "." + this.table1 + ".bare.wig";
		PrintWriter pwAs = new PrintWriter(new FileOutputStream(outfnAs));
		pwAs.printf("track type=wiggle_0 name=%sbare description=%sbare  color=" + colorA + " visibility=full " +
				" graphType=points autoScale=off alwaysZero=off maxHeightPixels=64:32:10 viewLimits=0:100\n", this.table1, this.table1);
		openFiles.add(pwAs);
	
		String outfnBs = this.outPrefix + "." + this.table2 + ".bare.wig";
		PrintWriter pwBs = new PrintWriter(new FileOutputStream(outfnBs));
		pwBs.printf("track type=wiggle_0 name=%sbare description=%sbare color=" + colorB + " visibility=full " +
				" graphType=points autoScale=off alwaysZero=off maxHeightPixels=64:32:10 viewLimits=0:100\n", this.table2, this.table2);
		openFiles.add(pwBs);
		
		
		
		String outfnA = this.outPrefix + "." + this.table1 + ".wig";
		PrintWriter pwA = new PrintWriter(new FileOutputStream(outfnA));
		pwA.printf("track type=wiggle_0 name=%s description=%s  color=" + colorA + " visibility=full " +
				" autoScale=off alwaysZero=off maxHeightPixels=64:32:10 viewLimits=0:100\n", this.table1, this.table1);
		openFiles.add(pwA);
	
		String outfnB = this.outPrefix + "." + this.table2 + ".wig";
		PrintWriter pwB = new PrintWriter(new FileOutputStream(outfnB));
		pwB.printf("track type=wiggle_0 name=%s description=%s color=" + colorB + " visibility=full " +
				" autoScale=off alwaysZero=off maxHeightPixels=64:32:10 viewLimits=0:100\n", this.table2, this.table2);
		openFiles.add(pwB);
		
		
		
		// HEMI
		String outfnHemiA = this.outPrefix + ".hemimeth." + this.table1 + ".wig";
		PrintWriter pwHemiA = new PrintWriter(new FileOutputStream(outfnHemiA));
		pwHemiA.printf("track type=wiggle_0 name=hemi-%s description=hemi-%s  color=204,153,102 visibility=full " +
				" autoScale=off alwaysZero=off maxHeightPixels=64:64:10 viewLimits=-50:50\n", this.table1, this.table1);
		openFiles.add(pwHemiA);
	
		String outfnHemiB = this.outPrefix + ".hemimeth." + this.table2 + ".wig";
		PrintWriter pwHemiB = new PrintWriter(new FileOutputStream(outfnHemiB));
		pwHemiB.printf("track type=wiggle_0 name=hemi-%s description=hemi-%s color=204,102,0 visibility=full " +
				" autoScale=off alwaysZero=off maxHeightPixels=64:64:10 viewLimits=-50:50\n", this.table2, this.table2);
		openFiles.add(pwHemiB);

		String outfnHemifakeA = this.outPrefix + ".hemifakemeth." + this.table1 + ".wig";
		PrintWriter pwHemifakeA = new PrintWriter(new FileOutputStream(outfnHemifakeA));
		pwHemifakeA.printf("track type=wiggle_0 name=hemifake-%s description=hemifake-%s  color=204,153,102 visibility=full " +
				" autoScale=off alwaysZero=off maxHeightPixels=64:64:10 viewLimits=-50:50\n", this.table1, this.table1);
		openFiles.add(pwHemifakeA);
	
		String outfnHemifakeB = this.outPrefix + ".hemifakemeth." + this.table2 + ".wig";
		PrintWriter pwHemifakeB = new PrintWriter(new FileOutputStream(outfnHemifakeB));
		pwHemifakeB.printf("track type=wiggle_0 name=hemifake-%s description=hemifake-%s color=204,102,0 visibility=full " +
				" autoScale=off alwaysZero=off maxHeightPixels=64:64:10 viewLimits=-50:50\n", this.table2, this.table2);
		openFiles.add(pwHemifakeB);

		
		String outfnDplus = this.outPrefix + ".deltaplus.wig";
		PrintWriter pwDplus = new PrintWriter(new FileOutputStream(outfnDplus));
		pwDplus.printf("track type=wiggle_0 name=%s description=%s color=255,0,0 visibility=full " + 
				" autoScale=off alwaysZero=off maxHeightPixels=64:32:10 viewLimits=0:60\n", "deltaplus", "deltaplus");
		openFiles.add(pwDplus);

		String outfnDminus = this.outPrefix + ".deltaMinus.wig";
		PrintWriter pwDminus = new PrintWriter(new FileOutputStream(outfnDminus));
		pwDminus.printf("track type=wiggle_0 name=%s description=%s color=0,255,0 visibility=full " + 
				" autoScale=off alwaysZero=off maxHeightPixels=64:32:10 viewLimits=-60:0\n", "deltaMinus", "deltaMinus");
		openFiles.add(pwDminus);

		
		String outCvg = this.outPrefix + ".cvgDelta.wig";
		PrintWriter pwCvg = new PrintWriter(new FileOutputStream(outCvg));
		pwCvg.printf("track type=wiggle_0 name=%s description=%s visibility=full " + 
				" autoScale=off alwaysZero=off maxHeightPixels=64:48:10 viewLimits=-2:2\n", "CoverageDelta", "CoverageDelta");
		openFiles.add(pwCvg);

		
		for (String chr : chrs)
		{

			if (!chrEndSet || (chrEnd<=0))
			{
				chrEnd = GoldAssembly.chromLengthStatic(chr, "hg18");
			}
			
			String rawspanline = String.format("variableStep chrom=%s\n",chr);
			pwAs.append(rawspanline);
			pwBs.append(rawspanline);

			//			
			// This is wildly ineffeficient, but i'm in a hurry right now.  Implement a
			// walker class in the future.
			int JUMP_INTERVAL = (int)Math.round((double)this.windSize/1.0);
			int MAXCPGS = 10000;
			List<String> tables = Arrays.asList(this.table1,this.table2);
			int mpCount = 0;
			int numCpgs = 0, flank=0, windS = 0, windE = 0;
			double[] abuffer = new double[MAXCPGS]; 
			double[] bbuffer = new double[MAXCPGS]; 
			double[] afwbuffer = new double[MAXCPGS]; 
			double[] bfwbuffer = new double[MAXCPGS]; 
			double[] arevbuffer = new double[MAXCPGS]; 
			double[] brevbuffer = new double[MAXCPGS]; 
			double[] afwfakebuffer = new double[MAXCPGS]; 
			double[] bfwfakebuffer = new double[MAXCPGS]; 
			double[] arevfakebuffer = new double[MAXCPGS]; 
			double[] brevfakebuffer = new double[MAXCPGS]; 
			double[] dbuffer = new double[MAXCPGS]; 
			double[] cvgabuffer = new double[MAXCPGS]; 
			double[] cvgbbuffer = new double[MAXCPGS]; 
			final double LOG2 = Math.log(2.0);

			boolean inSpan = false;
			boolean[] posSeen = new boolean[250000000];
			for (int mp = chrSt; mp <= chrEnd; mp += this.stepSize, mpCount++)
			{
				if ((mpCount % 1E2)==0) System.err.printf("Chr %s\tOn mp=%d (#%d)\n",chr, mp, mpCount);


				// We start at our window size, then go out until we have minCpgs.  We quit
				// if we hit maxWindStretch
				boolean done = false;
				for (int windCur = this.windSize; !done && (windCur <= this.maxWindStretch); windCur += JUMP_INTERVAL)
				{
					flank = windCur/2;
					windS = mp - flank;
					windE = mp + flank;

					params.clearRangeFilters();
					params.addRangeFilter(chr, windS, windE);
					CpgIteratorMultisample cpgit = new CpgIteratorMultisample(params, tables);
					numCpgs = cpgit.getCurNumRows();


					if (numCpgs >= this.minCpgs)
					{
						done = true;
						if ((mpCount % 1E2)==0) System.err.printf("\tFound a window windCur=%d (%d-%d): %d cpgs\n", windCur, windS, windE, numCpgs);

						// Start a new span sec if necessary
						if (!inSpan)
						{
							inSpan = true;
							String spanline = String.format("fixedStep chrom=%s start=%d step=%d span=%d\n",
									chr, mp, this.stepSize, (this.stepSize*2)-2);
							pwA.append(spanline);
							pwB.append(spanline);
							pwHemiA.append(spanline);
							pwHemiB.append(spanline);
							pwHemifakeA.append(spanline);
							pwHemifakeB.append(spanline);
							pwDplus.append(spanline);
							pwDminus.append(spanline);
							pwCvg.append(spanline);
						}


						// Now calculate the values!
						int numCpgsCounted = 0, numCpgsCountedFw = 0, numCpgsCountedRev = 0, numCpgsCountedFwfake = 0, numCpgsCountedRevfake = 0;
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

							double metha = cpgs[0].fracMeth(!this.noNonconvFilter);
							double methb = cpgs[1].fracMeth(!this.noNonconvFilter);
							if (Double.isNaN(metha) || Double.isNaN(methb))
							{
								System.err.printf("Why is meth value of CpG at %d NaN when i used a selective query??\n\t\t%s\n\t\t%s\n",
										cpgs[0].chromPos, cpgs[0].toStringExpanded(), cpgs[1].toStringExpanded());
							}
							else
							{
								// Finally we can do what we came to do .
								double delta = methb-metha;
								abuffer[numCpgsCounted] = metha;
								bbuffer[numCpgsCounted] = methb;
								dbuffer[numCpgsCounted] = delta;
								cvgabuffer[numCpgsCounted] = cpgs[0].totalReads+cpgs[0].totalReadsOpposite;
								cvgbbuffer[numCpgsCounted] = cpgs[1].totalReads+cpgs[1].totalReadsOpposite;
								numCpgsCounted++;

								// Strand specific
								if (cpgs[0].getStrand() == StrandedFeature.NEGATIVE)
								{
									arevbuffer[numCpgsCountedRev] = metha;
									brevbuffer[numCpgsCountedRev] = methb;
									numCpgsCountedRev++;
								}
								else
								{
									afwbuffer[numCpgsCountedFw] = metha;
									bfwbuffer[numCpgsCountedFw] = methb;
									numCpgsCountedFw++;
								}
								
								// Fake strand specific
								if (Math.random()>0.5)
								{
									arevfakebuffer[numCpgsCountedRevfake] = metha;
									brevfakebuffer[numCpgsCountedRevfake] = methb;
									numCpgsCountedRevfake++;
								}
								else
								{
									afwfakebuffer[numCpgsCountedFwfake] = metha;
									bfwfakebuffer[numCpgsCountedFwfake] = methb;
									numCpgsCountedFwfake++;
								}
								
								
								// And add the raw ones straight away
								int pos = cpgs[0].chromPos;
								if (this.bare && !posSeen[pos])
								{
									pwAs.printf("%d\t%d\n",pos,(int)Math.round(100.0*metha));
									pwBs.printf("%d\t%d\n",pos,(int)Math.round(100.0*methb));
									posSeen[pos] = true;
								}
							}
						} // foreach Cpgs

						// Write line
						int aval = (int)Math.round( 100.0 * MatUtils.nanMean(abuffer, 0, numCpgsCounted));
						pwA.println(aval);
						int bval = (int)Math.round( 100.0 * MatUtils.nanMean(bbuffer, 0, numCpgsCounted));
						pwB.println(bval);

						int ahemival = (int)Math.round( 100.0 * (MatUtils.nanMean(afwbuffer, 0, numCpgsCountedFw) -
								MatUtils.nanMean(arevbuffer, 0, numCpgsCountedRev)));
						pwHemiA.println(ahemival);
						
						int bhemival = (int)Math.round( 100.0 * (MatUtils.nanMean(bfwbuffer, 0, numCpgsCountedFw) -
								MatUtils.nanMean(brevbuffer, 0, numCpgsCountedRev)));
						pwHemiB.println(bhemival);

						int ahemifakeval = (int)Math.round( 100.0 * (MatUtils.nanMean(afwfakebuffer, 0, numCpgsCountedFwfake) -
								MatUtils.nanMean(arevfakebuffer, 0, numCpgsCountedRevfake)));
						pwHemifakeA.println(ahemifakeval);
						
						int bhemifakeval = (int)Math.round( 100.0 * (MatUtils.nanMean(bfwfakebuffer, 0, numCpgsCountedFwfake) -
								MatUtils.nanMean(brevfakebuffer, 0, numCpgsCountedRevfake)));
						pwHemifakeB.println(bhemifakeval);

						
						int dval = (int)Math.round( 100.0 * MatUtils.nanMean(dbuffer, 0, numCpgsCounted));
						pwDplus.println( (dval>0) ? dval : "0");
						pwDminus.println( (dval<0) ? dval : "0");

						// 0.9 is to correct for our actual disparity in coverage right now
						//					double cvgRatio = Math.log(MatUtils.nanMean(cvgbbuffer, 0, numCpgsCounted) / 
						//							(0.9 * MatUtils.nanMean(cvgabuffer, 0, numCpgsCounted)))/LOG2; 

						double cvg = MatUtils.nanMean(cvgbbuffer, 0, numCpgsCounted);
						pwCvg.printf("%.3f\n",cvg);

					} // (numCpgs>=minCpgs)

					cpgit.destroy();
				} // loop through window sizes

				if (!done)
				{
					System.err.printf("\tFound a window WITH INSUFFICIENT CPGs (%d-%d): %d cpgs\n", windS, windE, numCpgs);
					inSpan = false;
				}


			}
		}

		for (PrintWriter pw : openFiles)
		{
			pw.close();
		}

		
	} // Main

	public static double safeVal(double v)
	{
		return ((Double.isNaN(v) || Double.isInfinite(v)) ? 0.0 : v);
	}
	
}



//If memory gets to be a problem with whole genome loads, try this trick to walk the chromosome.		
//// Stupid JDBC tries to load entire chromosome into memory at once,
//// which is too much.
//final int STEP = (int)1E7;
//final int MAXCOORD = (int)2.8E8;
//for (int c = 0; c < MAXCOORD; c += STEP)
//{
//
//	Logger.getLogger(Logger.GLOBAL_LOGGER_NAME).severe("Starting chrom " + chrom + " from " + c+1 + " to " + c+STEP + "\n");
//
//	querier.clearRangeFilters();
//	querier.addRangeFilter(chrom, c+1, c+STEP);		

