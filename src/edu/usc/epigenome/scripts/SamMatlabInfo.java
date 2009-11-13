package edu.usc.epigenome.scripts;

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

import edu.usc.epigenome.genomeLibs.MiscUtils;
import edu.usc.epigenome.genomeLibs.PicardUtils;
import edu.usc.epigenome.genomeLibs.Counters.SNPByCycleCounter;


public class SamMatlabInfo {

	static final String USAGE = "SamMatlabInfo -summaryInfo -mismatchesByCycle -mismatchContexts -debug -minMapQ 30 MD XM ...";



	
	/**
	 * @param args
	 */

	@Option(name="-minMapQ",usage="minimum mapping quality (default 0)")
	protected int minMapQ = 0;
	@Option(name="-summaryInfo",usage=" An output which gives summary info for each alignment (default false)")
	protected boolean summaryInfo = false;
	@Option(name="-mismatchesByCycle",usage=" An alternative output which gives a table of mismatches, stratified by cycle (default false)")
	protected boolean mismatchesByCycle = false;
	@Option(name="-mismatchContexts",usage=" An alternative output which shows mismatches by context, one base before, one after (default false)")
	protected boolean mismatchContexts = false;
	@Option(name="-debug",usage=" Debugging statements (default false)")
	protected boolean debug = false;

	// receives other command line parameters than options
	@Argument
	private List<String> extraFields = new ArrayList<String>();

	

	
	/**
	 * @param args
	 */
	public static void main(String[] args)
	throws Exception
	{
		new SamMatlabInfo().doMain(args);
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
		
		
		
		final SAMFileReader inputSam = new SAMFileReader(System.in);
		inputSam.setValidationStringency(SAMFileReader.ValidationStringency.SILENT);
		int recCounter = 0;
		SNPByCycleCounter counter = new SNPByCycleCounter();
		
		record: for (final SAMRecord samRecord : inputSam) {
			// Convert read name to upper case.
			//samRecord.setReadName(samRecord.getReadName().toUpperCase());
						

			int mapQual = samRecord.getMappingQuality();
			boolean unmapped = samRecord.getReadUnmappedFlag();
			
			if (unmapped || (mapQual < minMapQ))
			{
				continue record;
			}
			
			recCounter++;
			if ((recCounter % 1E3)==0) System.err.println("On new record #" + recCounter); 

			

			try
			{

				for (String field : extraFields)  
				{
					Object val = samRecord.getAttribute(field);
					if (val == null) val = "nan";  // For MATLAB
					if (summaryInfo) System.out.print(val + ",");
				}
				
				if (summaryInfo) System.out.print(mapQual + ",");


				
				String md = (String)samRecord.getAttribute("MD");
				Cigar cigar = samRecord.getCigar();
				String seq = samRecord.getReadString().toUpperCase();
				String baseQual = samRecord.getBaseQualityString();
//				System.err.println("\tMD=" + md);
//				System.err.println("\tseq=" + seq);
				String ref = PicardUtils.refStrFromMd(seq, md, cigar).toUpperCase();

				// Revcomp if negative strand
				boolean revStrand = samRecord.getReadNegativeStrandFlag();
				if (revStrand)
				{
					seq = MiscUtils.revCompNucStr(seq);
					ref = MiscUtils.revCompNucStr(ref);
				}
				
				
				if (summaryInfo) System.out.print(seq + "\t");
				if (summaryInfo) System.out.print(baseQual + ",");
				if (summaryInfo) System.out.println();
				
				if (debug) System.err.printf("%s\t%d\t%s\t%s\n%s\n\n",seq, mapQual, md, samRecord.getCigarString(), ref);

//				if (seq.length() != ref.length())
//				{
//					System.err.println("SeqLen(" + seq.length() + ") != RefLen(" + ref.length() + ")");
//					System.err.println(seq + "\n" + ref);
//				}

				int seqLen = Math.min(seq.length(), ref.length());
				for (int i = 0; i < seqLen; i++)
				{
					char seqi = seq.charAt(i);
					char refi = ref.charAt(i);
					if (refi != '0')
					{
						if (debug && (seqi != refi)) System.err.printf("Mismatch %c->%c pos %d\n",refi,seqi,i);
						if (mismatchesByCycle) counter.increment(refi, seqi, i);
					}

				}
			}
			catch (Exception e)
			{
				System.err.println("-----------------------------------------");
				System.err.println("Couldn't handle seq #" + recCounter);
				System.err.println(samRecord.getReadName());
				e.printStackTrace(System.err);
				System.err.println("-----------------------------------------");

				System.exit(1);
			}

		} // record

		if (mismatchesByCycle)
		{
			System.out.print(counter.toString());
		}
		
		
		inputSam.close();

		
	}

	
	
//	static boolean isCytosine(int pos, String refStr, String seqStr)
//	{
//		char refC = refStr.charAt(pos);
//		char seqC = seqStr.charAt(pos);
//		
//		return ((refC == 'C') && ((seqC == 'C') || (seqC == 'T'))); 
//	}
//	
//	static boolean isCpg(int pos, String refStr, String seqStr)
//	{
//		if (pos >= (refStr.length()-1)) return false; // At the last character
//		
//		char refCnext = refStr.charAt(pos+1);
//		char seqCnext = seqStr.charAt(pos+1);
//		
//		return ( isCytosine(pos,refStr,seqStr) && (refCnext == 'G') && (seqCnext == 'G') );
//	}
//	
//	static boolean isConverted(int pos, String refStr, String seqStr)
//	{
//		char refC = refStr.charAt(pos);
//		char seqC = seqStr.charAt(pos);
//		
//		return ((refC == 'C') && (seqC == 'T'));
//	}
//
//	static int boolToInd(boolean bool)
//	{
//		return ( bool ? 0 : 1 );
//	}
//
//	static String headerLine()
//	{
//		return "\tpre\t\t\t\tpost" + "\n" + "\tcpg\t\tcph\t\tcpg\t\tcph" + "\n" + "convS\tc\tnc\tc\tnc\tc\tnc\tc\tnc"; 
//	}
//	
	
	
	

	
}
