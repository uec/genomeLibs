package edu.usc.epigenome.scripts;

import net.sf.samtools.*;
import net.sf.samtools.util.CloseableIterator;

import java.io.File;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.SortedMap;
import java.util.TreeMap;
import java.util.logging.Logger;

import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;

import edu.usc.epigenome.genomeLibs.MiscUtils;
import edu.usc.epigenome.genomeLibs.PicardUtils;
import edu.usc.epigenome.genomeLibs.MethylDb.Cpg;
import edu.usc.epigenome.genomeLibs.MethylDb.MethylDbUtils;


public class SamToMethyldbOffline {

	final private static String prefix = "methylCGsRich_";
	final private static String USAGE = "SamToMethyldbOffline [opts] sampleName file1.bam file2.bam ...";

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
	@Option(name="-outputCphs",usage=" Output CpH cytosines (can't use more than 1 input file)")
	protected boolean outputCphs = false;
	@Option(name="-minMapQ",usage="minimum mapping quality (default 0)")
	protected int minMapQ = 0;
	@Option(name="-minCphCoverage",usage="the minimum number of total reads to include a Cph (default 10)")
	protected int minCphCoverage = 10;
	@Option(name="-removeDups",usage="Removes reads starting and ending on the same position (false)")
	protected boolean removeDups = false;
	@Option(name="-allowBadMates",usage="Allows incorrectly mated paired ends (false)")
	protected boolean allowBadMates = false;
	@Option(name="-minCphFrac",usage="minimum methylation fraction to include a Cph")
	protected double minCphFrac = 0.2;
	@Option(name="-maxCoverageOutput",usage="maximum coverage for output array (default 100)")
	protected int maxCoverageOutput = 100;
//	@Option(name="-minOppStrandCoverage",usage="minimum coverage of opposite strand to confirm cytosine (default 5)")
//	protected int minOppStrandCoverage = 5;
//	@Option(name="-minNextBaseCoverage",usage="minimum coverage of next base to distinguish CpG/CpH (default 5)")
//	protected int minNextBaseCoverage = 5;
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
		new SamToMethyldbOffline().doMain(args);
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
		
		long totalUniqueMappedBasesRead = 0;
		short[] refCoordsCovered = new short[300000000]; // This really assumes a single chromosome
		
		// Iterate through chroms
		if (chrs.size()==0) chrs = MethylDbUtils.CHROMS;
		for (final String chr : chrs)
		{
			SortedMap<Integer,Cpg> cytosines = new TreeMap<Integer,Cpg>();
			PrintWriter outWriter = Cpg.outputChromToFile(cytosines, prefix, sampleName, chr, this.minCphCoverage, this.minCphFrac);

			
			for (final String fn : stringArgs)
			{
				File inputSamOrBamFile = new File(fn);

				final SAMFileReader inputSam = new SAMFileReader(inputSamOrBamFile);
				inputSam.setValidationStringency(SAMFileReader.ValidationStringency.SILENT);
				
				int querys = 10243;
				int querye = 10244;
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
						//System.err.printf("On base %d, purging everything before %d\n", lastBaseSeen,lastPurge);
						Cpg.outputCpgsToFile(outWriter, cytosines.headMap(new Integer(lastPurge)), prefix, sampleName, chr, this.minCphCoverage, this.minCphFrac);
						
						// Weird, if i just set cytosines to be the tailMap (as in 1, below) garbage collection doesn't actually clean up
						// the old part (just backed by the original data structure). So i actually have to copy it to a new map. (as in 2)
						//(1) cytosines = cytosines.tailMap(new Integer(lastPurge));
						cytosines = new TreeMap<Integer,Cpg>(cytosines.tailMap(new Integer(lastPurge))); // (2)
						
						lastPurge = lastBaseSeen;
					}
					
					SAMRecord samRecord = chrIt.next();
					//System.err.printf("New rec pos = %d\n",samRecord.getAlignmentStart());
					
					// Filter low qual
					int mapQual = samRecord.getMappingQuality();
					boolean unmapped = samRecord.getReadUnmappedFlag();
					if (unmapped || (mapQual < minMapQ))
					{
						// System.err.printf("Read unmapped or minQ<%d\n",minMapQ);
						continue record;
					}
					
					// Filter on uniqueness
					if (samRecord.getNotPrimaryAlignmentFlag())
					{
//						System.err.printf("Not primary alignment\n");
						continue record;
					}
					
					// Inverted dups, count only one end
					if (samRecord.getAlignmentStart() == samRecord.getMateAlignmentStart() && samRecord.getReadNegativeStrandFlag() == samRecord.getMateNegativeStrandFlag())
					{
						if (samRecord.getSecondOfPairFlag()) continue record;
						//System.err.printf("Inverted dup %d%s (%s)\n", samRecord.getAlignmentStart(), samRecord.getReadNegativeStrandFlag()?"-":"+", PicardUtils.getReadString(samRecord, true));
					}

					// If it's paired-end, filter on good mate unless otherwise specified
					if (samRecord.getReadPairedFlag()  && !allowBadMates && !samRecord.getProperPairFlag())
					{
//						System.err.printf("Not proper pair\n");
						continue record;
					}

					
					
					int alignmentS = samRecord.getAlignmentStart();
					if (this.removeDups)
					{
						if (samRecord.getDuplicateReadFlag())
						{
							//System.err.printf("Remove dups1\n");
							continue record;
						}
						
						// Compatibility with MAQ
						if (!samRecord.getReadPairedFlag())
						{
							if (lastBaseSeen == alignmentS) 
							{
								//System.err.printf("Remove dups2\n");
								continue record;
							}
						}
						
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
						
						
						// If we're in paired-end mode, and we're on the second end, we need to basically
						// reverse complement everything to get back to the cytosine strand.
						// The big problem with this is that 5' bisulfite conversion filter won't work right.
						boolean secondOfPair = samRecord.getReadPairedFlag() && getSecondOfPair(samRecord); // samRecord.getSecondOfPairFlag()
						if (secondOfPair)
						{
							negStrand = !negStrand;
							seq = MiscUtils.revCompNucStr(seq);
							ref = MiscUtils.revCompNucStr(ref);
						}
						
						boolean startingNeg = negStrand;
						if (secondOfPair) startingNeg = !startingNeg; 
						int	onRefCoord = (startingNeg) ? samRecord.getUnclippedEnd() : alignmentS;

						//System.err.printf("Seq (%d-%d,%s,%s) startRefCoord=%d: %s\n %s\n",samRecord.getAlignmentStart(), samRecord.getAlignmentEnd(), (negStrand?"-":"+"), (secondOfPair?"-":"+"), onRefCoord, seq,ref);

						
						if ((this.outputCphs) && (alignmentS < lastBaseSeen))
						{
							System.err.printf("BAM must be ordered in order to handle -outputCphs: %d<%d\n",alignmentS, lastBaseSeen);
							System.exit(1);
						}
						lastBaseSeen = alignmentS;
							
							

						int numConverted = 0;
						int convStart = Integer.MAX_VALUE;
						int seqLen = Math.min(seq.length(), ref.length());
						// If we're the second of the pair, we do this backwards so the conversion filter will work properly.
						int iStart = (secondOfPair) ? (seqLen-1) : 0;
						int iEnd = (secondOfPair) ? 0 : (seqLen - 1);
						int iInc = (secondOfPair) ? -1 : 1;
						
//						for (int i = 0; i < seqLen; i++)
						for (int i = iStart; i != iEnd; i += iInc)
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
								boolean iscpg = PicardUtils.isCpg(i,ref); // Rely only on reference cytosines
								boolean conv = PicardUtils.isConverted(i,ref,seq);
								
								//if (iscpg) System.err.printf("CpG at %d (%s%s,%s%s)\n", onRefCoord, refi, nextBaseRef, seqi, nextBaseSeq);
								

								if (conv && (this.useCpgsToFilter || !iscpg)) numConverted++;

								// If this is the first legal one , note it
								if ((convStart==Integer.MAX_VALUE) && (numConverted>=this.minConv) )
								{
									convStart = i;
								}


								boolean inNonconversionZone = (convStart==Integer.MAX_VALUE) || (secondOfPair&&(i>convStart)) || (!secondOfPair && (i<convStart));
								if (iscpg || this.outputCphs)
								{
									if (inNonconversionZone)
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
									
									this.incrementCpg(cpg, seqi, inNonconversionZone, nextBaseSeq);
								}
								
								// Count CpHs as a byproduct
								if (!iscpg)
								{
									//if (secondOfPair && inNonconversionZone) System.err.printf("secondOfPair & inNonconversion: %s:%d\n",chr,onRefCoord);
									
									numCphTotalNoFilt++;
									if (conv) numCphConvertedNoFilt++;
									
									if (!inNonconversionZone)
									{
										numCphTotalWithFilt++;
										if (conv) numCphConvertedWithFilt++;
									}
								}



							} // IsCytosine

//							// Since we are doing CpHs, we need opp strand info for all guanines
//							if (PicardUtils.isGuanine(i, ref))
							boolean oppositeCpg = PicardUtils.isOppositeCpg(i,ref);
							if (oppositeCpg)
							{
								// Look for cpg on opposite strand
								// Be careful, the "nextBaseRef" is now on the opposite strand!!
								char nextBaseRefRev = PicardUtils.nextBaseRef(i, ref, true);
								Cpg cpg = findOrCreateCpg(cytosines, onRefCoord, !negStrand, nextBaseRefRev);
								
								this.incrementOppositeCpg(cpg, seqi);
							}


							// Increment genomic coord position. Careful, we go backwards if we're
							// on the second of a pair.
							if (refi == '-')
							{
								// It's a deletion in reference, don't advance
							}
							else
							{
								boolean backwardsRef = negStrand;
								if (secondOfPair) backwardsRef = !backwardsRef; // We do not do this! We walk in the same direction as the original bisulfite strand! // 
								int inc = (backwardsRef) ? -1 : 1;
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

			// And output the chrom
			if (cytosines != null) 	Cpg.outputCpgsToFile(outWriter, cytosines, prefix, sampleName, chr, this.minCphCoverage, this.minCphFrac);
			outWriter.close();

		}
		double frac = (double)cpgsFilteredOutCounter/((double)cpgsUsedCounter+(double)cpgsFilteredOutCounter);
		System.err.printf("Lost %f%% due to non-converion filter\n%d CpGs filtered for non-conversion, %d CpGs used (MinConv=%d,UseCpgs=%s)\n",
				frac*100.0, cpgsFilteredOutCounter, cpgsUsedCounter, this.minConv, String.valueOf(this.useCpgsToFilter));
		System.err.printf("Found %d reads total\n", recCounter);
		System.err.printf("Cph conversion rate: before filter=%f, after filter=%f\n",
				100.0 * ((double)numCphConvertedNoFilt/(double)numCphTotalNoFilt),
				100.0 * ((double)numCphConvertedWithFilt/(double)numCphTotalWithFilt));
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
	
	
	
	protected static boolean getSecondOfPair(SAMRecord read) {
		return read.getSecondOfPairFlag();
//		boolean secondOfPair = false;
//		String readName = read.getReadName();
//	    final String END1_SUFFIX = String.format("%c1", '/');
//	    final String END2_SUFFIX = String.format("%c2", '/');
//		if (readName.endsWith(END1_SUFFIX))
//		{
//			secondOfPair = false;
//		}
//		else if (readName.endsWith(END2_SUFFIX))
//		{
//			secondOfPair = true;   			
//		}
//		else
//		{
//			System.err.println("Got a read that doesn't end with /1 or /2: " + readName + ".  Can't tell which end it is.");
//			System.exit(1);
//		}	
//		
//		return secondOfPair;
	}
}
