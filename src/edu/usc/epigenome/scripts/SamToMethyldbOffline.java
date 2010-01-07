package edu.usc.epigenome.scripts;

import net.sf.samtools.*;
import net.sf.samtools.util.CloseableIterator;

import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;

import edu.usc.epigenome.genomeLibs.PicardUtils;
import edu.usc.epigenome.genomeLibs.MethylDb.Cpg;
import edu.usc.epigenome.genomeLibs.MethylDb.MethylDbUtils;


public class SamToMethyldbOffline {

	final private static String prefix = "methylCGsRich_";
	final private static String USAGE = "SamToMethyldbOffling [opts] sampleName file1.bam file2.bam ...";

	
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
		int usedCounter = 0;
		int filteredOutCounter = 0;

		// Iterate through chroms
		if (chrs.size()==0) chrs = MethylDbUtils.CHROMS;
		for (final String chr : chrs)
		{
			Map<Integer,Cpg> cpgs = new TreeMap<Integer,Cpg>();

			for (final String fn : stringArgs)
			{
				File inputSamOrBamFile = new File(fn);

				final SAMFileReader inputSam = new SAMFileReader(inputSamOrBamFile);
				inputSam.setValidationStringency(SAMFileReader.ValidationStringency.SILENT);
				
				CloseableIterator<SAMRecord> chrIt = inputSam.query(chr, 0, 0, false);
				

				record: while (chrIt.hasNext())
				{
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
					if ((recCounter % 1E5)==0) System.err.println("On new record #" + recCounter); 


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
						int	onRefCoord = (negStrand) ? samRecord.getUnclippedEnd() : samRecord.getAlignmentStart();


						int numConverted = 0;
						int convStart = Integer.MAX_VALUE;
						int seqLen = Math.min(seq.length(), ref.length());
						for (int i = 0; i < seqLen; i++)
						{
							char refi = ref.charAt(i);
							char seqi = seq.charAt(i);

							if ((seqi != '0') && PicardUtils.isCytosine(i,ref))
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

									Cpg cpg = findOrCreateCpg(cpgs, onRefCoord, negStrand);
									this.incrementCpg(cpg, seqi, i<convStart);
								}



							} // IsCytosine

							boolean oppositeCpg = PicardUtils.isOppositeCpg(i,ref);
							if (oppositeCpg)
							{
								// Look for cpg on opposite strand
								Cpg cpg = findOrCreateCpg(cpgs, onRefCoord, !negStrand);
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
						chrIt.close();
						System.exit(1);
					}

				} // record

				chrIt.close();
				inputSam.close();

			}

			// And output the chrom
			if (cpgs != null) Cpg.outputChromToFile(cpgs, prefix, sampleName, chr);

		}
		double frac = (double)filteredOutCounter/((double)usedCounter+(double)filteredOutCounter);
		System.err.printf("Lost %f%% due to non-converion filter\n%d CpGs filtered for non-conversion, %d CpGs used (MinConv=%d,UseCpgs=%s)\n",
				frac*100.0, filteredOutCounter, usedCounter, this.minConv, String.valueOf(this.useCpgsToFilter));
		System.err.printf("Found %d reads total\n", recCounter);
	}

	protected static Cpg findOrCreateCpg(Map<Integer,Cpg> cpgs, int onRefCoord, boolean negStrand)
	{
		Cpg cpg = cpgs.get(new Integer(onRefCoord));
		if (cpg == null)
		{
			cpg = new Cpg(onRefCoord,negStrand);
			cpgs.put(new Integer(onRefCoord), cpg);
		}
		return cpg;
	}
	
	
	protected void incrementCpg(Cpg cpg, char seqChar, boolean nonconvFilter) 
	throws Exception
	{
		int totalReads = 0, cReads = 0, tReads = 0, cReadsNonconvFilt = 0, agReads = 0;
		
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

		cpg.totalReads += totalReads;
		cpg.cReads += cReads;
		cpg.tReads += tReads;
		cpg.cReadsNonconversionFilt += cReadsNonconvFilt;
		cpg.agReads += agReads;
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
