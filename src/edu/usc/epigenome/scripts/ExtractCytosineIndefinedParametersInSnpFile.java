package edu.usc.epigenome.scripts;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;

import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.util.CloseableIterator;

import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;

import edu.usc.epigenome.genomeLibs.PicardUtils;


public class ExtractCytosineIndefinedParametersInSnpFile {
	/**
	 * @param args
	 */
		// it is not correct to use for SNP array site for CpH and CpG...
		//final private static String prefix = "methylCGsRich_ASM_AllSnp_";
		final private static String USAGE = "ExtractCytosineIndefinedParametersInSnpFile [opts] snpFileName bamFile ...";

		final private static int PURGE_INTERVAL = 20000; // We purge our stored Cpgs once we get this many bases past them.
		//final private static int ALLELE_GA_NUMBER = 1; //for each reads, when there are more than 1 GA position is different from reference sequence, we define it belongs to another allele. 
		
		/**
		 * @param args
		 */
		
		@Option(name="-minMapQ",usage="minimum mapping quality (default 30)")
		protected int minMapQ = 30;
		@Option(name="-minReadCov",usage="minimum read coverage (default 1)")
		protected int minReadCov = 1;
		@Option(name="-minCTReadCov",usage="minimum CT read coverage (default 3)")
		protected int minCTReadCov = 3;
		@Option(name="-minBaseQual",usage="minimum Base quality (default 0)")
		protected int minBaseQual = 0;
		@Option(name="-minOppositeGCount",usage="minimum G Allele Count (default 3)")
		protected static double minOppositeGCount = 3;
		
		@Option(name="-minOppositeACount",usage="minimum A Allele Count (default 1)")
		protected static double minOppositeACount = 1;
		
		@Option(name="-minOppositeAFreq",usage="minimum BAllele Frequency (default 0.1)")
		protected static double minOppositeAFreq = 0.10;
		@Option(name="-minCTContent",usage="minCT read in the C site (default 0.9)")
		protected static double minCTContent = 0.90;
		@Option(name="-minGuanContent",usage="minGuan read in the CpG sites' G site (default 0.9)")
		protected static double minGuanContent = 0.90;
		@Option(name="-minHContent",usage="minGuan read in the CpH sites' H site (default 0.9)")
		protected static double minHContent = 0.90;
		//@Option(name="-cpg",usage=" seperately consider CpG and CpH sites (default false)")
		//protected boolean cpg = false;
		//@Option(name="-female",usage="the sample is female or male (default false for male)")
		//protected boolean female = false;
		@Option(name="-debug",usage=" Debugging statements (default false)")
		protected boolean debug = false;
		// receives other command line parameters than options
		@Argument
		private List<String> stringArgs = new ArrayList<String>();

		// when readCov<10, a>=1 will be filt out, if readCov>=10, a>10% will be filt out
		
		/**
		 * @param args
		 */
		public static void main(String[] args)
		throws Exception
		{
			ExtractCytosineIndefinedParametersInSnpFile test = new ExtractCytosineIndefinedParametersInSnpFile();
			test.doMain(args);
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


			String snpFileName = stringArgs.remove(0);
			String bamFile = stringArgs.remove(0);
			BufferedReader br = new BufferedReader(new FileReader(snpFileName));
			String line;
			
			File inputSamOrBamFile = new File(bamFile);
			SAMFileReader inputSam = new SAMFileReader(inputSamOrBamFile);
			inputSam.setValidationStringency(SAMFileReader.ValidationStringency.SILENT);
			
			String cytosineFileName = snpFileName + ".minReadCov_" + minReadCov + ".minCTReadCov_" + minCTReadCov + ".cytosinecaller.all.txt";
			PrintWriter cytosineWriter = new PrintWriter(new File(cytosineFileName));
			int hom = 0;
			int het = 0;
			

			while( (line = br.readLine()) != null){
				String[] tmpArray = line.split("\t");
			//	if(!tmpArray[5].equalsIgnoreCase("C") && !tmpArray[6].equalsIgnoreCase("C") && !tmpArray[5].equalsIgnoreCase("G") && !tmpArray[6].equalsIgnoreCase("G"))
				//	continue;
				
				
				String chr = "chr" + tmpArray[1];
				
				
				
				Integer snpPosition =(int) Math.round(Double.parseDouble(tmpArray[2]));
				
				boolean cytosineFlag = false;
				
				
				
				int totalNum = 0;
				
				int totalNumPos = 0;

				int totalNumNeg = 0;
	
				int totalNumC = 0;
				int totalNumT = 0;

			
				
				
				CloseableIterator<SAMRecord> chrIt = inputSam.queryOverlapping(chr, snpPosition, snpPosition);
				record: while (chrIt.hasNext())
				{	
					SAMRecord samRecord = chrIt.next();
					int mapQual = samRecord.getMappingQuality();
					byte[] baseQual = samRecord.getBaseQualities();
					boolean unmapped = samRecord.getReadUnmappedFlag();
					if (unmapped || (mapQual < minMapQ))
					{
						continue record;
					}
					String seq = PicardUtils.getReadString(samRecord, true);
					try{
						
						String ref = PicardUtils.refStr(samRecord, true);
						boolean negStrand = samRecord.getReadNegativeStrandFlag();
						int alignmentS = samRecord.getAlignmentStart();
						
						int	onRefCoord = (negStrand) ? samRecord.getUnclippedEnd() : alignmentS; 
						int i = Math.abs(snpPosition - onRefCoord);
						if (seq.length() != ref.length())
						{
							System.err.println("SeqLen(" + seq.length() + ") != RefLen(" + ref.length() + ")");
							System.err.println(seq + "\n" + ref);
						}
						int seqLen = Math.min(seq.length(), ref.length());
						byte baseQS = (negStrand) ? baseQual[seqLen-1-i] : baseQual[i];
						if( baseQS < minBaseQual )
						{
							continue record;
						}
						if (ref.charAt(i) == '-')
						{
							break record;
						}
						
						totalNum++;
						if(negStrand){
							totalNumNeg++;
					
						}
						else{
							totalNumPos++;
							
						}
	
						if(PicardUtils.isCytosine(i, seq, true)){
							cytosineFlag=true;
							if(PicardUtils.isCytosine(i, seq, false)){
								totalNumC++;
							}
							else{
								totalNumT++;
							}
							
								
							
						}
						
						
					}
					catch (Exception e)
					{
						System.err.println("-----------------------------------------");
						System.err.println("Couldn't handle seq #");
						System.err.println(seq);
						e.printStackTrace(System.err);
						System.err.println("-----------------------------------------");
//						chrIt.close();
//						System.exit(1);
					}
					
				}
				chrIt.close();
				if(totalNum < minReadCov || (totalNumC + totalNumT) < minCTReadCov)
					continue;
				cytosineWriter.println(line);
				if(tmpArray[5].equalsIgnoreCase(tmpArray[6]))
					hom++;
				else
					het++;
				
			}	
			cytosineWriter.close();
			br.close();
			inputSam.close();
			System.out.println("Het: " + het + "\tHom: " + hom);
	}

}
