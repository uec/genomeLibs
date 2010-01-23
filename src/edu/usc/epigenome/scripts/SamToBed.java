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

import edu.usc.epigenome.genomeLibs.PicardUtils;
import edu.usc.epigenome.genomeLibs.MethylDb.Cpg;
import edu.usc.epigenome.genomeLibs.MethylDb.MethylDbUtils;


public class SamToBed {

	final private static String prefix = "methylCGsRich_";
	final private static String USAGE = "SamToBed [opts] sampleName file1.bam file2.bam ... > out.bed";

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
//	@Option(name="-numCycles",usage="Number of cycles to track")
//	protected int numCycles = 100;
//	@Option(name="-outputReads",usage=" Outputs one line per read (default false)")
//	protected boolean outputReads = false;
	@Option(name="-useCpgsToFilter",usage=" Use CpGs and CpHs to filter if true, otherwise just CpHs (default false)")
	protected boolean useCpgsToFilter = false;
	@Option(name="-minMapQ",usage="minimum mapping quality (default 0)")
	protected int minMapQ = 0;
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
		new SamToBed().doMain(args);
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
		
		System.out.printf("track name=\"%s\" description=\"%s\" useScore=0 itemRgb=On visibility=4\n",
				sampleName, sampleName);

		
		int recCounter = 0;
		int usedCounter = 0;
		int filteredOutCounter = 0;

		// Iterate through chroms
		if (chrs.size()==0) chrs = MethylDbUtils.CHROMS;
		for (final String chr : chrs)
		{
			for (final String fn : stringArgs)
			{
				File inputSamOrBamFile = new File(fn);

				final SAMFileReader inputSam = new SAMFileReader(inputSamOrBamFile);
				inputSam.setValidationStringency(SAMFileReader.ValidationStringency.SILENT);
				
				CloseableIterator<SAMRecord> chrIt = inputSam.query(chr, 0, 0, false);
				
				record: while (chrIt.hasNext())
				{

					SAMRecord samRecord = chrIt.next();
					
					// Use a consensus CpG
					Cpg consensusCpg = null;

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
						System.err.printf("On new record #%d\n",recCounter);
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
						int	onRefCoord = (negStrand) ? samRecord.getUnclippedEnd() : alignmentS; 

							

						int numConverted = 0;
						int convStart = Integer.MAX_VALUE;
						int seqLen = Math.min(seq.length(), ref.length());
						for (int i = 0; i < seqLen; i++)
						{
							char refi = ref.charAt(i);
							char seqi = seq.charAt(i);
							char nextBaseRef = PicardUtils.nextBaseRef(i, ref);
							char nextBaseSeq = PicardUtils.nextBaseSeq(i, seq);

							//if ((i < (seqLen-1)) && PicardUtils.isCytosine(i,ref)) // The last one is too tricky to deal with since we don't know context
							if ((i < (seqLen-1)) && PicardUtils.isCytosine(i,ref,false) && PicardUtils.isCytosine(i, seq,true)) // The last one is too tricky to deal with since we don't know context
							{
								boolean iscpg = PicardUtils.isCpg(i,ref);
								boolean conv = PicardUtils.isConverted(i,ref,seq);
																

								if (conv && (this.useCpgsToFilter || !iscpg)) numConverted++;

								// If this is the first legal one , note it
								if ((convStart==Integer.MAX_VALUE) && (numConverted>=this.minConv) )
								{
									convStart = i;
								}


								if (iscpg)
								{
									if (i<convStart)
									{
										// In the non-conversion filter zone
										filteredOutCounter++;
										//System.err.printf("Rec %d\tpos=%d\n",recCounter,i);
									}
									else
									{
										// Past the non-conversion filter, use it
										usedCounter++;
									}

									if (consensusCpg == null)
									{
										consensusCpg = createCpg(onRefCoord, negStrand, nextBaseRef);
									}
									
									// See if we can fix the context for this CpG
									if (consensusCpg.getNextBaseRef() == '0') consensusCpg.setNextBaseRef(nextBaseRef);
									
									this.incrementCpg(consensusCpg, seqi, i<convStart, nextBaseSeq);

								}



							} // IsCytosine


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

					// Output bed line
					String line = toBedLine(samRecord, consensusCpg);
					if (line != null) System.out.println(line);
					
					
				} // record

				chrIt.close();
				inputSam.close();

			} // BAM file


		}
		double frac = (double)filteredOutCounter/((double)usedCounter+(double)filteredOutCounter);
		System.err.printf("Lost %f%% due to non-converion filter\n%d CpGs filtered for non-conversion, %d CpGs used (MinConv=%d,UseCpgs=%s)\n",
				frac*100.0, filteredOutCounter, usedCounter, this.minConv, String.valueOf(this.useCpgsToFilter));
		System.err.printf("Found %d reads total\n", recCounter);
	}

	/**
	 * @param rec
	 * @param consensusCpg
	 * @return return null if not enough cpgs
	 */
	public String toBedLine(SAMRecord rec, Cpg consensusCpg)
	{
		
		String strand = rec.getReadNegativeStrandFlag() ? "-" : "+";
		int s = rec.getAlignmentStart();
		int e = rec.getAlignmentEnd();
		String chrom = rec.getReferenceName();
		
		double meth = (consensusCpg == null) ? Double.NaN : consensusCpg.fracMeth(true);

		int score;
		String color = methToColor(meth);
		if (Double.isNaN(meth) || Double.isInfinite(meth))
		{
			score = 0;
		}
		else if (meth > 0.5)
		{
			score = (int)Math.round(100.0 * ((meth-0.5)*2.0));
		}
		else
		{
			score = (int)Math.round(100.0 * ((0.5-meth)*2.0));
		}

		String out = String.format("%s\t%d\t%d\t%s\t%d\t%s\t%d\t%d\t%s\t", 
				chrom,
				s,
				e,
				String.format("r%d%s", s,strand),
				score,
				strand,  // strand information gets in the way of display
				s,
				e,
				color
				);
		return out;
	}

	protected static String methToColor(double meth)
	{
		String color="255,255,204"; // 

		if (Double.isNaN(meth))
		{
		}
		else if (meth > 0.5)
		{
			int dec = (int)Math.floor(10.0 * ((meth-0.5)*2.0));
			// Red
//			switch (dec)
//			{
//			case 0: color = "49,49,49"; break;
//			case 1: color = "55,44,44"; break;
//			case 2: color = "61,38,38"; break;
//			case 3: color = "66,33,33"; break;
//			case 4: color = "72,27,27"; break;
//			case 5: color = "78,22,22"; break;
//			case 6: color = "83,16,16"; break;
//			case 7: color = "89,11,11"; break;
//			case 8: color = "94,5,5"; break;
//			case 9: case 10: color = "100,0,0"; break;
//			default: System.err.println("Got illegal color decile: " + dec); System.exit(1); break;
//			}
			switch (dec)
			{
			case 0: color = "153,153,153"; break;
			case 1: case 2: case 3: color = "255,153,153"; break;
			case 4: case 5: case 6: case 7: color = "255,102,102"; break;
			case 8: case 9: case 10: color = "255,0,0"; break;
			default: System.err.println("Got illegal color decile: " + dec); System.exit(1); break;
			}
		}
		else
		{
			int dec = (int)Math.floor(10.0 * ((0.5-meth)*2.0));
			// Green
//			switch (dec)
//			{
//			case 0: color = "49,49,49"; break;
//			case 1: color = "44,55,44"; break;
//			case 2: color = "38,61,38"; break;
//			case 3: color = "33,66,33"; break;
//			case 4: color = "27,72,27"; break;
//			case 5: color = "22,78,22"; break;
//			case 6: color = "16,83,16"; break;
//			case 7: color = "11,89,11"; break;
//			case 8: color = "5,94,5"; break;
//			case 9: case 10: color = "0,100,0"; break;
//			default: System.err.println("Got illegal color decile: " + dec); System.exit(1); break;
//			}
			switch (dec)
			{
			case 0: color = "153,153,153"; break;
			case 1: case 2: case 3: color = "153,255,153"; break;
			case 4: case 5: case 6: case 7: color = "102,255,102"; break;
			case 8: case 9: case 10: color = "0,255,0"; break;
			default: System.err.println("Got illegal color decile: " + dec); System.exit(1); break;
			}
		}
	
		return color;
	}
	
	
	protected static Cpg createCpg(int onRefCoord, boolean negStrand, char nextBaseRef)
	{
		Cpg cpg = new Cpg(onRefCoord,negStrand);
		cpg.setNextBaseRef(nextBaseRef);

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
