package edu.usc.epigenome.testScripts;

import net.sf.samtools.*;
import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;

import edu.usc.epigenome.genomeLibs.PicardUtils;


public class SamFiveprimeConversion {

	static final String USAGE = "SamFiveprimeConversion -minConv 1 input.sam";

	static final int BUFFERLEN = 1000;
	static StringBuilder sb = new StringBuilder(BUFFERLEN);
	
	static Pattern pairPat = Pattern.compile("([0-9]*)([ACTGNactgn]?)");

	
	/**
	 * @param args
	 */

	@Option(name="-minConv",usage="minimum number of converted cytosines required")
	protected int minConv = 1;
	@Option(name="-numCycles",usage="Number of cycles to track")
	protected int numCycles = 100;
	@Option(name="-outputReads",usage=" Outputs one line per read (default false)")
	protected boolean outputReads = false;
	@Option(name="-useCpgsToFilter",usage=" Use CpGs and CpHs to filter if true, otherwise just CpHs (default false)")
	protected boolean useCpgsToFilter = false;

	// receives other command line parameters than options
	@Argument
	private List<String> stringArgs = new ArrayList<String>();

	
	public static void clearStrBuf()
	{
		sb.delete(0, BUFFERLEN);
	}
	
	/**
	 * @param args
	 */
	public static void main(String[] args)
	throws Exception
	{
		new SamFiveprimeConversion().doMain(args);
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
			if (stringArgs.size() != 1) throw new CmdLineException("No input file specified");
		}
		catch (CmdLineException e)
		{
			System.err.println(e.getMessage());
			System.err.println(USAGE);
			// print the list of available options
			parser.printUsage(System.err);
			System.err.println();
			return;
		}
		
		
		if (outputReads) System.out.println(headerLine());
		
		File inputSamOrBamFile = new File(stringArgs.get(0));

		CytosineCounter globalCounter = (new SamFiveprimeConversion()).new CytosineCounter();
		ByCycleCounter cpgCycleCounter = (new SamFiveprimeConversion()).new ByCycleCounter(numCycles);
		ByCycleCounter cphCycleCounter = (new SamFiveprimeConversion()).new ByCycleCounter(numCycles);
		

		final SAMFileReader inputSam = new SAMFileReader(inputSamOrBamFile);
		inputSam.setValidationStringency(SAMFileReader.ValidationStringency.SILENT);
		int recCounter = 0;
		record: for (final SAMRecord samRecord : inputSam) {
			// Convert read name to upper case.
			//samRecord.setReadName(samRecord.getReadName().toUpperCase());
						
			String seq = PicardUtils.getReadString(samRecord, true);

			recCounter++;
			if ((recCounter % 1E5)==0) System.err.println("On new record #" + recCounter); 
			
			try
			{

				String ref = PicardUtils.refStr(samRecord, true);
				//			System.err.println("\tref=" + ref);

				if (seq.length() != ref.length())
				{
					System.err.println("SeqLen(" + seq.length() + ") != RefLen(" + ref.length() + ")");
					System.err.println(seq + "\n" + ref);
				}

				CytosineCounter localCounter = (new SamFiveprimeConversion()).new CytosineCounter();

				int numConverted = 0;
				int convStart = Integer.MAX_VALUE;
				int seqLen = Math.min(seq.length(), ref.length());
				for (int i = 0; i < seqLen; i++)
				{
					if (isCytosine(i,ref,seq))
					{
						boolean cpg = isCpg(i,ref,seq);
						boolean conv = isConverted(i,ref,seq);

						if (conv && (useCpgsToFilter || !cpg)) numConverted++;

						// Don't count the first one, it throws off the stats
						if ((convStart==Integer.MAX_VALUE) && (numConverted>=minConv) )
						{
							convStart = i;
						}
						else
						{
							localCounter.increment(i<convStart, cpg, conv);
							globalCounter.increment(i<convStart, cpg, conv);

							ByCycleCounter cycleCounter = (cpg) ? cpgCycleCounter : cphCycleCounter;
							cycleCounter.increment(i<convStart, conv, i);
						}
						
						//if (cpg && (i<convStart)) System.err.printf("Rec %d\tpos=%d\n",recCounter,i);

					}
				}

				if (outputReads) System.out.println(convStart + "\t" + localCounter.toString() + "\t" + seq + "\t" + ref);
			}
			catch (Exception e)
			{
				System.err.println("-----------------------------------------");
				System.err.println("Couldn't handle seq #" + recCounter);
				System.err.println(seq);
				e.printStackTrace(System.err);
				System.err.println("-----------------------------------------");
			}

		} // record

		
		
		inputSam.close();
		System.out.println("\n" + headerLine() + "\n\t" + globalCounter.toString());
		
		System.out.println("\nCpG cycle counter");
		System.out.println(cpgCycleCounter.toString());

		System.out.println("\nCpH cycle counter");
		System.out.println(cphCycleCounter.toString());
		
	}

	
	
	static boolean isCytosine(int pos, String refStr, String seqStr)
	{
		char refC = refStr.charAt(pos);
		char seqC = seqStr.charAt(pos);
		
		return ((refC == 'C') && ((seqC == 'C') || (seqC == 'T'))); 
	}
	
	static boolean isCpg(int pos, String refStr, String seqStr)
	{
		if (pos >= (refStr.length()-1)) return false; // At the last character
		
		char refCnext = refStr.charAt(pos+1);
		char seqCnext = seqStr.charAt(pos+1);
		
		return ( isCytosine(pos,refStr,seqStr) && (refCnext == 'G') && (seqCnext == 'G') );
	}
	
	static boolean isConverted(int pos, String refStr, String seqStr)
	{
		char refC = refStr.charAt(pos);
		char seqC = seqStr.charAt(pos);
		
		return ((refC == 'C') && (seqC == 'T'));
	}

	static int boolToInd(boolean bool)
	{
		return ( bool ? 0 : 1 );
	}

	static String headerLine()
	{
		return "\tpre\t\t\t\tpost" + "\n" + "\tcpg\t\tcph\t\tcpg\t\tcph" + "\n" + "convS\tc\tnc\tc\tnc\tc\tnc\tc\tnc"; 
	}
	
	
	public class ByCycleCounter
	{
		// ind 1 = Pre or post.
		// ind 2 = conv/non-converted
		// ind 3 = cycle
		int[][][] counter; 
		int numCycles = 0;
		
		
		/**
		 * @param numCycles
		 */
		public ByCycleCounter(int numCycles) {
			this.numCycles = numCycles;
			counter = new int[2][2][numCycles];
		}

		void increment(boolean pre, boolean converted, int cycle)
		{
			if (cycle < numCycles)
			{
				counter[boolToInd(pre)][boolToInd(converted)][cycle]++;
			}
		}

		public String toString()
		{
			// Do it with a stringBuf so it's faster
			clearStrBuf();
			
			for (int i = 0; i < 2; i++)
			{
				for (int j = 0; j < 2; j++)
				{
					sb.append("byCycle line: " + i + ", " + j);
					for (int k = 0; k < numCycles; k++)
					{
						sb.append(counter[i][j][k] + ",");
					}
					sb.append("\n");
				}
			}			
			
			return sb.toString();
		}
	
	}
	
	public class CytosineCounter
	{
		
		// ind 1 = Pre or post.
		// ind 2 = CpG or CpH 
		// ind 3 = converted/non-converted
		int[][][] counter = new int[2][2][2];
		
		void increment(boolean pre, boolean cpg, boolean converted)
		{
			//System.err.println("Pre = " + pre);
			counter[boolToInd(pre)][boolToInd(cpg)][boolToInd(converted)]++;
		}
		
		public String toString()
		{
			// Do it with a stringBuf so it's faster
			clearStrBuf();
			
			for (int i = 0; i < 2; i++)
			{
				for (int j = 0; j < 2; j++)
				{
					for (int k = 0; k < 2; k++)
					{
						sb.append(counter[i][j][k] + "\t");
					}
				}
			}
			
			return sb.toString();
		}
	
	}
	
	
	
	
	/*
	 * Roll this back into Picard
	 * 
	 */
	static String refStrFromMd(String seq, String md)
	throws Exception
	{
		clearStrBuf();

		boolean done = false; 
		Matcher match = pairPat.matcher(md);
		int curSeqPos = 0;
		while (!done)
		{
			if (match.find())
			{
				int num = -1;
				char nuc = 0;
				String mg;
				if ((mg = match.group(1)).length() > 0)
				{
					num = Integer.parseInt(mg);
					for (int i = 0; i < num; i++)
					{
						sb.append(seq.charAt(curSeqPos++));
					}
				}
				if ((mg = match.group(2)).length() > 0)
				{
					nuc = mg.charAt(0);
					sb.append(nuc);
					curSeqPos++;
				}
				
				done = match.hitEnd();
//				System.err.println("\tnum=" + num + "\tnuc=" +nuc);
				
			}
			else
			{
				throw new Exception("Illegal MD pattern: " + md);
			}

		}
		
		return sb.toString();
	}
	
	
}
