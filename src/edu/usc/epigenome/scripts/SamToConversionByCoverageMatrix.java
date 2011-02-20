package edu.usc.epigenome.scripts;

import net.sf.samtools.*;
import net.sf.samtools.util.CloseableIterator;

import java.io.File;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.SortedMap;
import java.util.TreeMap;
import java.util.logging.Logger;

import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;

import edu.usc.epigenome.genomeLibs.ListUtils;
import edu.usc.epigenome.genomeLibs.MatUtils;
import edu.usc.epigenome.genomeLibs.PicardUtils;
import edu.usc.epigenome.genomeLibs.MethylDb.Cpg;
import edu.usc.epigenome.genomeLibs.MethylDb.CytosineStats;
import edu.usc.epigenome.genomeLibs.MethylDb.MethylDbUtils;


public class SamToConversionByCoverageMatrix {

	final private static String prefix = "conversionByCoverage_";
	final private static String USAGE = "SamToConversionByCoverageMatrix [opts] sampleName file1.bam file2.bam ...";

//	final private static int PURGE_INTERVAL = 100000; // We purge our stored Cpgs once we get this many bases past them.
	final private static int PURGE_INTERVAL = 20000; // We purge our stored Cpgs once we get this many bases past them.
	
	
	/**
	 * object vars
	 */

	
	/**
	 * @param args
	 */
	@Option(name="-chrom",multiValued=true,usage="One or more chroms, eg. --chrom chr1 --chrom chr5")
	protected List<String> chrs = new ArrayList<String>(25);
	@Option(name="-minConv",usage="minimum number of converted cytosines required")
	protected int minConv = 1;
	@Option(name="-useCpgsToFilter",usage=" Use CpGs and CpHs to filter if true, otherwise just CpHs (default false)")
	protected boolean useCpgsToFilter = false;
	@Option(name="-outputCphs",usage=" Output CpH cytosines (default false)")
	protected boolean outputCphs = false;
	@Option(name="-minMapQ",usage="minimum mapping quality (default 0)")
	protected int minMapQ = 0;
	@Option(name="-maxCoverageOutput",usage="maximum coverage for output array (default 100)")
	protected int maxCoverageOutput = 100;
	@Option(name="-minOppStrandCoverage",usage="minimum coverage of opposite strand to confirm cytosine (default 5)")
	protected int minOppStrandCoverage = 5;
	@Option(name="-minNextBaseCoverage",usage="minimum coverage of next base to distinguish CpG/CpH (default 5)")
	protected int minNextBaseCoverage = 5;
	@Option(name="-debug",usage=" Debugging statements (default false)")
	protected boolean debug = false;
	@Option(name="-maxOppStrandAfrac",usage="As on the opposite strand are evidence for mutation or SNP. " +
	"This sets a maximum number of observed As on the opposite strand (default 0.1)")
	protected double maxOppStrandAfrac = 0.101;
	@Option(name="-maxNextNonGfrac",usage="If the base following the C has more than this ratio of non-G bases, we don't count it. (default 0.1)")
	protected double maxNextNonGfrac = 0.101;
	
	// receives other command line parameters than options
	@Argument
	private List<String> stringArgs = new ArrayList<String>();

	
	
	/**
	 * @param args
	 */
	public static void main(String[] args)
	throws Exception
	{
		new SamToConversionByCoverageMatrix().doMain(args);
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
			if (stringArgs.size() < 2) throw new CmdLineException(USAGE);
			
			if (this.outputCphs && (stringArgs.size()>2)) throw new CmdLineException("Can't use -outputCphs with multiple input bam files");
		}
		catch (CmdLineException e)
		{
			System.err.println(e.getMessage());
			// print the list of available options
			parser.printUsage(System.err);
			System.err.println();
			return;
		}

		String sampleName = stringArgs.remove(0);

		int recCounter = 0;
		int cpgsUsedCounter = 0;
		int cpgsFilteredOutCounter = 0;

		// Cph counters
		int numCphConvertedWithFilt = 0;
		int numCphTotalWithFilt = 0;
		int numCphConvertedNoFilt = 0;
		int numCphTotalNoFilt = 0;
		
		int totalUniqueMappedBasesRead = 0;
		short[] refCoordsCovered = new short[300000000]; // This really assumes a single chromosome
		CytosineStats cytStats = new CytosineStats(this.maxOppStrandAfrac, this.maxNextNonGfrac);
		
		// Iterate through chroms
		if (chrs.size()==0) chrs = MethylDbUtils.CHROMS;
		for (final String chr : chrs)
		{
			// Output map
			int[][] cphCounts = new int[maxCoverageOutput][maxCoverageOutput]; 
			int[][] cpgCounts = new int[maxCoverageOutput][maxCoverageOutput]; 

			// Current cytosines
			SortedMap<Integer,Cpg> cytosines = new TreeMap<Integer,Cpg>();
			
			for (final String fn : stringArgs)
			{
				File inputSamOrBamFile = new File(fn);

				final SAMFileReader inputSam = new SAMFileReader(inputSamOrBamFile);
				inputSam.setValidationStringency(SAMFileReader.ValidationStringency.SILENT);
				
				int querys = 7000000;
				int querye = 29000000;
				querys = 0;
				querye = 0;
				CloseableIterator<SAMRecord> chrIt = inputSam.query(chr, querys, querye, false);
				
				// We can only purge if we are the only input file and we are sorted.
				boolean canPurge = ((inputSam.hasIndex()) && (stringArgs.size() == 1));
				Logger.getLogger(Logger.GLOBAL_LOGGER_NAME).info(String.format("Able to purge at interval %d? %s\n", 
						PURGE_INTERVAL, (canPurge) ? "Yes" : "Purging not available - either unsorted BAM or multiple BAMs for the same chrom"));
				int lastBaseSeen = 0;
				int lastPurge = 0;
				record: while (chrIt.hasNext())
				{
					if (canPurge && (lastBaseSeen > (lastPurge+PURGE_INTERVAL)))
					{
						// First, increment counters
						if (cytosines != null) incrementCytosineCounters(cytosines.headMap(new Integer(lastPurge)), cpgCounts, cphCounts,cytStats);
						//System.err.printf("Purging %d cs from %d-%d\n", cytosines.headMap(new Integer(lastPurge)).size(), lastPurge, lastBaseSeen);

						
						// Weird, if i just set cytosines to be the tailMap (as in 1, below) garbage collection doesn't actually clean up
						// the old part (just backed by the original data structure). So i actually have to copy it to a new map. (as in 2)
						//(1) cytosines = cytosines.tailMap(new Integer(lastPurge));
						cytosines = new TreeMap<Integer,Cpg>(cytosines.tailMap(new Integer(lastPurge))); // (2)
						
						lastPurge = lastBaseSeen;
					}
					
					SAMRecord samRecord = chrIt.next();

					// Filter low qual
					int mapQual = samRecord.getMappingQuality();
					boolean unmapped = samRecord.getReadUnmappedFlag();
					if (unmapped || (mapQual < minMapQ))
					{
						continue record;
					}


					String seq = PicardUtils.getReadString(samRecord, true);

					recCounter++;
					if ((recCounter % 1E5)==0)
					{
						System.err.printf("On new record #%d, lastBase=%d, purged tree size:%d\n",recCounter,lastBaseSeen, cytosines.size());
						if (canPurge) System.gc();
					}


					try
					{
						String ref = PicardUtils.refStr(samRecord, true);

						if (seq.length() != ref.length())
						{
							System.err.println("SeqLen(" + seq.length() + ") != RefLen(" + ref.length() + ")");
							System.err.println(seq + "\n" + ref);
						}
						//System.err.println(seq + "\n" + ref);

						boolean negStrand = samRecord.getReadNegativeStrandFlag();
						int alignmentS = samRecord.getAlignmentStart();
						int	onRefCoord = (negStrand) ? samRecord.getUnclippedEnd() : alignmentS; // This will be incremented as we go along in the read

						if ((this.outputCphs) && (alignmentS < lastBaseSeen))
						{
							System.err.printf("BAM must be ordered in order to handle -outputCphs: %d<%d\n",alignmentS, lastBaseSeen);
							System.exit(1);
						}
						lastBaseSeen = alignmentS;
							
							

						int numConverted = 0;
						int convStart = Integer.MAX_VALUE;
						int seqLen = Math.min(seq.length(), ref.length());
						for (int i = 0; i < seqLen; i++)
						{
							totalUniqueMappedBasesRead++;
							refCoordsCovered[onRefCoord]++;
							
							char refi = ref.charAt(i);
							char seqi = seq.charAt(i);
							char nextBaseRef = PicardUtils.nextBaseRef(i, ref);
							char nextBaseSeq = PicardUtils.nextBaseSeq(i, seq);

							// We only look at cytosines in the reference
							if ((i < (seqLen-1)) && PicardUtils.isCytosine(i,ref,false) && PicardUtils.isCytosine(i, seq,true)) // The last one is too tricky to deal with since we don't know context
							{
								boolean iscpg = PicardUtils.isCpg(i,ref,seq);
								boolean conv = PicardUtils.isConverted(i,ref,seq);
								
								

								if (conv && (this.useCpgsToFilter || !iscpg)) numConverted++;

								// If this is the first legal one , note it
								if ((convStart==Integer.MAX_VALUE) && (numConverted>=this.minConv) )
								{
									convStart = i;
								}


								if (iscpg || this.outputCphs)
								{
									if (i<convStart)
									{
										// In the non-conversion filter zone
										if (iscpg) cpgsFilteredOutCounter++;
										//System.err.printf("Rec %d\tpos=%d\n",recCounter,i);
									}
									else
									{
										// Past the non-conversion filter, use it
										if (iscpg) cpgsUsedCounter++;
									}

									// onRefCoord is incremented below
									Cpg cpg = findOrCreateCpg(cytosines, onRefCoord, negStrand, nextBaseRef);
									// See if we can fix the context for this CpG
									if (cpg.getNextBaseRef() == '0') cpg.setNextBaseRef(nextBaseRef);
									
									this.incrementCpg(cpg, seqi, i<convStart, nextBaseSeq);
								}
								
								// Count CpHs as a byproduct
								if (!iscpg)
								{
									numCphTotalNoFilt++;
									if (conv) numCphConvertedNoFilt++;
									
									if (i>=convStart)
									{
										numCphTotalWithFilt++;
										if (conv) numCphConvertedWithFilt++;
									}
								}



							} // IsCytosine

//							boolean oppositeCpg = PicardUtils.isOppositeCpg(i,ref);
//							if (oppositeCpg)
							// Since we are doing CpHs, we need opp strand info for all guanines
							if (PicardUtils.isGuanine(i, ref))
							{
								// Look for cpg on opposite strand
								// Be careful, the "nextBaseRef" is now on the opposite strand!!
								char nextBaseRefRev = PicardUtils.nextBaseRef(i, ref, true);
								Cpg cpg = findOrCreateCpg(cytosines, onRefCoord, !negStrand, nextBaseRefRev);
								
								this.incrementOppositeCpg(cpg, seqi);
							}


							// Increment genomic coord position
							if (refi == '-')
							{
								// It's a deletion in reference, don't advance
							}
							else
							{
								int inc = (negStrand) ? -1 : 1;
								onRefCoord += inc;
							}

						} // i (pos within read)
					}
					catch (Exception e)
					{
						System.err.println("-----------------------------------------");
						System.err.println("Couldn't handle seq #" + recCounter);
						System.err.println(seq);
						e.printStackTrace(System.err);
						System.err.println("-----------------------------------------");
//						chrIt.close();
//						System.exit(1);
					}

				} // record

				chrIt.close();
				inputSam.close();

			}

			// Take care of last cytosines
			if (cytosines != null) incrementCytosineCounters(cytosines, cpgCounts, cphCounts,cytStats);

			// And output the chrom
			String fn = String.format("%s%s-%s-minOppCvg%d-minNextCvg%d-maxOppA%.3f-maxNextNonG%.3f-%s.csv", 
					prefix,sampleName,"CpG",this.minOppStrandCoverage, this.minNextBaseCoverage, this.maxOppStrandAfrac, this.maxNextNonGfrac,chr);
			PrintWriter outWriter = new PrintWriter(new File(fn));
			MatUtils.matlabCsv(outWriter, MatUtils.intMatToDouble(cpgCounts));
			outWriter.close();

			if (this.outputCphs)
			{
				fn = String.format("%s%s-%s-minOppCvg%d-minNextCvg%d-maxOppA%.3f-maxNextNonG%.3f-%s.csv", 
						prefix,sampleName,"CpH",this.minOppStrandCoverage, this.minNextBaseCoverage, this.maxOppStrandAfrac, this.maxNextNonGfrac,chr);
				outWriter = new PrintWriter(new File(fn));
				MatUtils.matlabCsv(outWriter, MatUtils.intMatToDouble(cphCounts));
				outWriter.close();
			}
		}
		

		// Print some stats
		String chrSec = (chrs.size()==1) ? ("-" + chrs.get(0)) : "";
//		String fn = String.format("bamBisulfiteStats-%s-minMapQ%d-maxOppA%.3f-maxNextNonG%.3f%s.csv", 
//				sampleName, this.minMapQ, this.maxOppStrandAfrac, this.maxNextNonGfrac, chrSec); 
		String fn = String.format("bamBisulfiteStats-%s%s.csv", sampleName, chrSec); 
		PrintWriter outWriter = new PrintWriter(new File(fn));

		
//		double frac = (double)cpgsFilteredOutCounter/((double)cpgsUsedCounter+(double)cpgsFilteredOutCounter);
//		outWriter.printf("Lost %f%% due to non-converion filter\n%d CpGs filtered for non-conversion, %d CpGs used (MinConv=%d,UseCpgs=%s)\n",
//				frac*100.0, cpgsFilteredOutCounter, cpgsUsedCounter, this.minConv, String.valueOf(this.useCpgsToFilter));
//		outWriter.printf("Found %d reads total\n", recCounter);
//		outWriter.printf("Cph conversion rate: before filter=%f, after filter=%f\n",
//				100.0 * ((double)numCphConvertedNoFilt/(double)numCphTotalNoFilt),
//				100.0 * ((double)numCphConvertedWithFilt/(double)numCphTotalWithFilt));
//		
		
		
		List<String> headers = new ArrayList<String>();
		List<String> vals = new ArrayList<String>();

		headers.add("Chrom");
		vals.add(String.format("%s", chrs.get(0)));

		headers.add("minMapQ"); 
		vals.add(String.format("%d", this.minMapQ));

		headers.add("maxOppStrandAfrac"); 
		vals.add(String.format("%.3f", this.maxOppStrandAfrac));

		headers.add("maxNextNonGfrac"); 
		vals.add(String.format("%.3f", this.maxNextNonGfrac));

		headers.add("minFivePrimeConverted"); 
		vals.add(String.format("%d", this.minConv));

		headers.add("totalReads");
		vals.add(String.format("%d", recCounter));

		headers.add("totalUniqueMappedBasesRead");
		vals.add(String.format("%d", totalUniqueMappedBasesRead));
		
		if (chrs.size()==1) // Otherwise these are not valid
		{
			int totalRefCoordsCovered = 0;
			int totalBases = 0;
			for (int i = 0; i < refCoordsCovered.length; i++)
			{
				short b = refCoordsCovered[i];
				if (b>0) totalRefCoordsCovered++;
				totalBases += b;
			}
			headers.add("totalRefCoordsCovered");
			vals.add(String.format("%d", totalRefCoordsCovered));

			headers.add("totalBases");
			vals.add(String.format("%d", totalBases));
		}
		
		headers.add("cpgsFivePrimeFiltered");
		vals.add(String.format("%d", cpgsFilteredOutCounter));

		headers.add("cpgsNotFivePrimeFiltered");
		vals.add(String.format("%d", cpgsUsedCounter));

		headers.add("cphConversionPreFilter");
		vals.add(String.format("%.3f%%", 100.0 * ((double)numCphConvertedNoFilt/(double)numCphTotalNoFilt)));
		
		
		headers.addAll(cytStats.headerStrings());
		vals.addAll(cytStats.summaryStrings());



		
		
		outWriter.println(ListUtils.excelLine(headers));
		outWriter.println(ListUtils.excelLine(vals));
		
		outWriter.close();

	}


	private void incrementCytosineCounters(SortedMap<Integer, Cpg> cytosines,
			int[][] cpgCounts, int[][] cphCounts, CytosineStats cytStats) 
	{
		incrementCytosineCounter(cytosines, cpgCounts, false);
		if (this.outputCphs) incrementCytosineCounter(cytosines, cphCounts, true);
		
		Iterator<Cpg> cytosineIt = cytosines.values().iterator();
		CPG: while (cytosineIt.hasNext())
		{
			Cpg cytosine = cytosineIt.next();
			cytStats.streamCytosine(cytosine);
		}
	}

	/**
	 * @param cytosines
	 * @param cpgCounts
	 * @param countCphs If true, we count only CpHs, if false we count only CpGs
	 */
	private void incrementCytosineCounter(SortedMap<Integer, Cpg> cytosines,
			int[][] cpgCounts, boolean countCphs) {

//		System.err.printf("Incrementing counter (cph=%s), total pre = %d\n",countCphs,MatUtils.sumAll(cpgCounts));
		Iterator<Cpg> cytosineIt = cytosines.values().iterator();
		int cytosineCount = 0;
		int goodCount = 0;
		int cpgCount = 0;
		int cphCount = 0;
		CPG: while (cytosineIt.hasNext())
		{
			Cpg cytosine = cytosineIt.next();
			cytosineCount++;

			//System.err.printf("\tCytosine %s\n",cytosine.toStringExpanded());
			
			
			// Only count it if it looks like a true cytosine.  Use a minimum of 5 reads for both cytosine
			// and CpG determination
			if ((cytosine.totalReadsCorT(true) >=1) && (cytosine.fracReadsCorT()>(1-this.maxOppStrandAfrac)) && 
					(cytosine.totalReadsOpposite>=this.minOppStrandCoverage) && (cytosine.fracOppositeA() <= this.maxOppStrandAfrac))
			{
				//System.err.printf("\t\tCytosine opp cvg>=5 %s\n",cytosine.toStringExpanded());

				if (cytosine.nextBaseTotalReads >= this.minNextBaseCoverage)
				{
					//System.err.printf("\t\tCytosine next base cvg>=5 %s\n",cytosine.toStringExpanded());

					boolean trueCpg = (cytosine.fracNextBaseG() >= (1.0-this.maxNextNonGfrac));
					boolean trueCph = (cytosine.fracNextBaseG() < this.maxNextNonGfrac);
					
					if (trueCpg) cpgCount++;
					if (trueCph) cphCount++;

					// If it's the right kind, count it
					if ((trueCpg && !countCphs) || (trueCph && countCphs))
					{
						if (cytosine.totalReadsCorT(true) < this.maxCoverageOutput)
						{
							cpgCounts[cytosine.totalReadsCorT(true)][cytosine.totalReadsC(true)]++;
							goodCount++;
						}
					}
				}
			}
		}
//		System.err.printf("\tIncremented counter (cph=%s), good cytosines=%d/%d, cpg=%d, cph=%d, total post = %d\n\n",
//				countCphs,goodCount,cytosineCount,cpgCount, cphCount, MatUtils.sumAll(cpgCounts));

		
	}

	protected static Cpg findOrCreateCpg(Map<Integer,Cpg> cpgs, int onRefCoord, boolean negStrand, char nextBaseRef)
	{
		Cpg cpg = cpgs.get(new Integer(onRefCoord));
		if (cpg == null)
		{
			cpg = new Cpg(onRefCoord,negStrand);
			cpg.setNextBaseRef(nextBaseRef);
			cpgs.put(new Integer(onRefCoord), cpg);
		}
		return cpg;
	}
	
	
	protected void incrementCpg(Cpg cpg, char seqChar, boolean nonconvFilter, char nextBaseSeq) 
	throws Exception
	{
		int totalReads = 0, cReads = 0, tReads = 0, cReadsNonconvFilt = 0, agReads = 0, nextBaseGreads = 0, nextBaseTotalReads = 0;
		
		switch (seqChar)
		{
		case 'N':
		case '0':
			break;
		case 'A':
		case 'G':
			agReads = 1;
			totalReads = 1;
			break;
		case 'T':
			tReads = 1;
			totalReads = 1;
			break;
		case 'C':
			if (nonconvFilter) cReadsNonconvFilt = 1; else cReads = 1;
			totalReads = 1;
			break;
		default:
			throw new Exception("Can't recognize seq char: " + seqChar);
		}
	
		// And the next base
		switch (nextBaseSeq)
		{
		case 'N':
		case '0':
			break;
		case 'G':
			nextBaseTotalReads = 1;
			nextBaseGreads = 1;
			break;
		case 'T':
		case 'A':
		case 'C':
			nextBaseTotalReads = 1;
			break;
		default:
			throw new Exception("Can't recognize seq char: " + seqChar);
		}
	
		cpg.totalReads += totalReads;
		cpg.cReads += cReads;
		cpg.tReads += tReads;
		cpg.cReadsNonconversionFilt += cReadsNonconvFilt;
		cpg.agReads += agReads;
		
		cpg.nextBaseGreads += nextBaseGreads;
		cpg.nextBaseTotalReads += nextBaseTotalReads;
	}
	
	protected void incrementOppositeCpg(Cpg cpg, char seqChar) 
	throws Exception
	{
		int aReadsOpposite = 0, totalReadsOpposite = 0;
		
		
		switch (seqChar)
		{
		case 'N':
		case '0':
			break;
		case 'A':
			aReadsOpposite = 1;
			totalReadsOpposite = 1;
			break;
		case 'G':
		case 'T':
		case 'C':
			totalReadsOpposite = 1;
			break;
		default:
			throw new Exception("Can't recognize seq char: " + seqChar);
		}

		cpg.aReadsOpposite += aReadsOpposite;
		cpg.totalReadsOpposite += totalReadsOpposite;
	}
	

	
}
