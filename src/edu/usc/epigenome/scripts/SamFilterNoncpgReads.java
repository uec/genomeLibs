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
import org.usckeck.genome.ListUtils;

import edu.usc.epigenome.genomeLibs.MiscUtils;
import edu.usc.epigenome.genomeLibs.PicardUtils;
import edu.usc.epigenome.genomeLibs.Counters.SNPByCycleCounter;


public class SamFilterNoncpgReads {

	static final String USAGE = "SamFilterNoncpgReads -dryRun -debug -minMapQ 30 file1.sam file2.bam ...";


	/**
	 * @param args
	 */

	@Option(name="-minMapQ",usage="minimum mapping quality (default 0)")
	protected int minMapQ = 0;
	@Option(name="-dryRun",usage=" Don't write output file")
	protected boolean dryRun = false;
	@Option(name="-debug",usage=" Debugging statements (default false)")
	protected boolean debug = false;

	// receives other command line parameters than options
	@Argument
	private List<String> inputFns = new ArrayList<String>();




	/**
	 * @param args
	 */
	public static void main(String[] args)
	throws Exception
	{
		new SamFilterNoncpgReads().doMain(args);
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
			if (inputFns.size() < 1) throw new CmdLineException("No input files specified.");
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

		// Start work
		fn: for (String inputFn : inputFns)
		{

			final SAMFileReader inputSam = new SAMFileReader(new File(inputFn));
			inputSam.setValidationStringency(SAMFileReader.ValidationStringency.SILENT);

			SAMFileWriter outputSam = null;
			if (!dryRun)
			{
				String outputFn = inputFn;
				outputFn = outputFn.replaceAll(".[sb]am", "");
				System.err.println("Out file name= " + outputFn);
				outputFn += ".onlyCpgReads.bam";
				outputSam = new SAMFileWriterFactory().makeSAMOrBAMWriter(inputSam.getFileHeader(),
						true, new File(outputFn));
			}
			
			
			int nTotal = 0, nMapped = 0, nCpg = 0;
			record: for (final SAMRecord samRecord : inputSam) {
				// Convert read name to upper case.
				//samRecord.setReadName(samRecord.getReadName().toUpperCase());

				nTotal++;
				if ((nTotal % 1E3)==0) System.err.println("On new record #" + nTotal); 


				// Filter low qual
				int mapQual = samRecord.getMappingQuality();
				boolean unmapped = samRecord.getReadUnmappedFlag();
				if (unmapped || (mapQual < minMapQ))
				{
					continue record;
				}
				nMapped++;


				try
				{

					boolean containsCpg = PicardUtils.readContainsCpg(samRecord);		
					if (containsCpg)
					{
						if (!dryRun) outputSam.addAlignment(samRecord);
						nCpg++;
					}
					


				}
				catch (Exception e)
				{
					System.err.println("-----------------------------------------");
					System.err.println("Couldn't handle seq #" + nTotal);
					System.err.println(samRecord.getReadName());
					e.printStackTrace(System.err);
					System.err.println("-----------------------------------------");

					System.exit(1);
				}

			} // record

			System.err.printf("File %s: nSeqs=%d, nMapped=%d, nContainCpg=%d\n", inputFn, nTotal, nMapped, nCpg);

			if (!dryRun) outputSam.close();
			inputSam.close();


		} // Input FN

	}





}
