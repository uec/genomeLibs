package edu.usc.epigenome.scripts;

import net.sf.samtools.*;
import net.sf.samtools.util.CloseableIterator;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;

import edu.usc.epigenome.genomeLibs.PicardUtils;
import edu.usc.epigenome.genomeLibs.MethylDb.Cpg;
import edu.usc.epigenome.genomeLibs.MethylDb.MethylDbUtils;


public class SamToMethylreaddbOffline {

	final private static String USAGE = "SamToMethylreaddbOffline [opts] file1.bam file2.bam ... > out.bed";

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
		new SamToMethylreaddbOffline().doMain(args);
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
			if (stringArgs.size() < 1) throw new CmdLineException(USAGE);
			
		}
		catch (CmdLineException e)
		{
			System.err.println(e.getMessage());
			// print the list of available options
			parser.printUsage(System.err);
			System.err.println();
			return;
		}

		
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
					
					// Use a consensus CpG to represent the whole read
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
					String line = toMysqlRow(samRecord, consensusCpg);
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
	public String toMysqlRow(SAMRecord rec, Cpg consensusCpg)
	{
		
		String strand = rec.getReadNegativeStrandFlag() ? "-" : "+";
		int s = rec.getAlignmentStart();
		int e = rec.getAlignmentEnd();
		String chrom = rec.getReferenceName();
		
		int nCpg = (consensusCpg == null) ? 0 : (consensusCpg.cReads + consensusCpg.tReads);
		int nCpgMeth = (consensusCpg == null) ? 0 : consensusCpg.cReads;
		int nCpgMajor = (consensusCpg == null) ? 0 : Math.max(consensusCpg.cReads, consensusCpg.tReads);
		
		// Only return lines with CpGs
		String out = null;
		if (nCpg > 0)
		{
			out = String.format("%d\t%d\t%s\t%d\t%d\t%d",
				s,
				e,
				strand,
				nCpg,
				nCpgMeth,
				nCpgMajor);
		}

		return out;
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
