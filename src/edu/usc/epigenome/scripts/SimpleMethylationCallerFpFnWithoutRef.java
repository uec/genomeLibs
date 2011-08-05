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

public class SimpleMethylationCallerFpFnWithoutRef {

	
		// it is not correct to use for SNP array site for CpH and CpG...
		//final private static String prefix = "methylCGsRich_ASM_AllSnp_";
		final private static String USAGE = "SimpleMethylationCallerFpFnWithoutRef [opts] snpFileName bamFile ...(without reference C)";

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
		@Option(name="-minOppositeGCount",usage="minimum G Allele frequency in the non bisulfite conversion strand (default 0.90)")
		protected static double minOppositeGCount = 0.90;
		
	//	@Option(name="-minOppositeACount",usage="minimum A Allele Count (default 1)")
	//	protected static double minOppositeACount = 1;
		
	//	@Option(name="-minOppositeAFreq",usage="minimum BAllele Frequency (default 0.1)")
	//	protected static double minOppositeAFreq = 0.10;
		@Option(name="-minCTContent",usage="minCT read in the bisulfite conversion strand (default 0.9)")
		protected static double minCTContent = 0.90;
	//	@Option(name="-minGuanContent",usage="minGuan read in the CpG sites' G site (default 0.9)")
	//	protected static double minGuanContent = 0.90;
	//	@Option(name="-minHContent",usage="minGuan read in the CpH sites' H site (default 0.9)")
	//	protected static double minHContent = 0.90;
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
			SimpleMethylationCallerFpFnWithoutRef test = new SimpleMethylationCallerFpFnWithoutRef();
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
			//int seqHeterLoci = 0;
			//int seqHomoLoci = 0;
			int snpHeterLociC = 0;
			int snpHomoLociC = 0;
			int snpHeterLociCafterCTreadcovFilter = 0;
			int snpHomoLociCafterCTreadcovFilter = 0;

			int truePC = 0;
			int falsePC = 0;
			int trueNC = 0;
			int falseNC = 0;
			
			
			
			int recCounter = 0;
			//String fn = sampleName + "Merge_chr1.NODUPS.sorted.calmd.NODUPS.bam";
			File inputSamOrBamFile = new File(bamFile);
			SAMFileReader inputSam = new SAMFileReader(inputSamOrBamFile);
			inputSam.setValidationStringency(SAMFileReader.ValidationStringency.SILENT);
			
			//String cytosineFileName = snpFileName + ".minReadCov_" + minReadCov + ".minCTReadCov_" + minCTReadCov + ".minGOpposit_" + minOppositeGCount + ".cytosinecaller.test.txt";
			//PrintWriter cytosineWriter = new PrintWriter(new File(cytosineFileName));
			
			
			
			while( (line = br.readLine()) != null){
				String[] tmpArray = line.split("\t");
				
				
				String chr = "chr" + tmpArray[1];
				recCounter++;
				if ((recCounter % 10000)==0)
				{
					System.err.printf("On new record #%d\n",recCounter);
					inputSam.close();
					
					System.gc();
					inputSam = new SAMFileReader(inputSamOrBamFile);
					inputSam.setValidationStringency(SAMFileReader.ValidationStringency.SILENT);
					
				}
				
				Integer snpPosition =(int) Math.round(Double.parseDouble(tmpArray[2]));
				
				boolean snpHomo = tmpArray[3].equalsIgnoreCase(tmpArray[4]);
				
				if((tmpArray[5].equalsIgnoreCase("A") && tmpArray[6].equalsIgnoreCase("A")) || (tmpArray[5].equalsIgnoreCase("T") && tmpArray[6].equalsIgnoreCase("T")) ||  (tmpArray[5].equalsIgnoreCase("A") && tmpArray[6].equalsIgnoreCase("T")) || (tmpArray[5].equalsIgnoreCase("T") && tmpArray[6].equalsIgnoreCase("A"))){
					
				}
				else{
					if(!snpHomo){
						snpHeterLociC++;
					}
					else{
						snpHomoLociC++;
					}
				}
				
				boolean cytosineFlag = false;
			//	boolean cytosineNegFlag = false;
				boolean cytosinePassFlag = false;
				
				
				int totalNum = 0;
				
				int totalNumPos = 0;

				int totalNumNeg = 0;
	
				int totalNumC = 0;
				int totalNumT = 0;

				int totalNumCTPos = 0;
				int totalNumCTNeg = 0;
				
				int totalNumGPos = 0;
				int totalNumGNeg = 0;
				
				
				
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
						//int readsStart = (negStrand) ? 0 - samRecord.getUnclippedEnd() : alignmentS;
						//int readsEnd = (negStrand) ? 0 - alignmentS : samRecord.getUnclippedEnd();
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
					//		refNegRecord = ref.charAt(i);
						//	if(seq.charAt(i) == 'A')
						//		totalNumCTNeg++;
							if(seq.charAt(i) == 'G')
								totalNumGNeg++;
						//	if(i < seqLen-1)
							//	totalNumNegnext++;
						}
						else{
							totalNumPos++;
							//refPosRecord = ref.charAt(i);
						//	if(seq.charAt(i) == 'A')
						//		totalNumCTPos++;
							if(seq.charAt(i) == 'G')
								totalNumGPos++;
						//	if(i < seqLen-1)
							//	totalNumPosnext++;
						}
	
						if(PicardUtils.isCytosine(i, seq, true)){
							cytosineFlag=true;
							if(PicardUtils.isCytosine(i, seq, false)){
								totalNumC++;
							}
							else{
								totalNumT++;
							}
							if(negStrand){
								//cytosineNegFlag = true;
								totalNumCTNeg++;
						
							}
							else{
								//cytosineNegFlag = false;
								totalNumCTPos++;
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
				if(totalNum < minReadCov && (totalNumC + totalNumT) < minCTReadCov)
					continue;
				
				if((tmpArray[5].equalsIgnoreCase("A") && tmpArray[6].equalsIgnoreCase("A")) || (tmpArray[5].equalsIgnoreCase("T") && tmpArray[6].equalsIgnoreCase("T")) ||  (tmpArray[5].equalsIgnoreCase("A") && tmpArray[6].equalsIgnoreCase("T")) || (tmpArray[5].equalsIgnoreCase("T") && tmpArray[6].equalsIgnoreCase("A"))){
					
				}
				else{
					if(!snpHomo){
						snpHeterLociCafterCTreadcovFilter++;
					}
					else{
						snpHomoLociCafterCTreadcovFilter++;
					}
				}
				
				if(cytosineFlag){
					if(((double)totalNumGPos/(double)totalNumPos > minOppositeGCount && (double)totalNumCTNeg/(double)totalNumNeg > minCTContent) || ((double)totalNumGNeg/(double)totalNumNeg > minOppositeGCount && (double)totalNumCTPos/(double)totalNumPos > minCTContent)){
						cytosinePassFlag = true;
					}
					
				}
				
				/*
				if(cytosineFlag){
					if(cytosineNegFlag){
						if(totalNumPos <= 0 && ((double)(totalNumC + totalNumT)/(double)totalNumNeg) >= minCTContent){
							cytosinePassFlag = true;
							
							
						}
						else if(totalNumPos > 0 && totalNumPos < 10 && ((double)(totalNumC + totalNumT)/(double)totalNumNeg) >= minCTContent){
							if(totalNumAPos < minOppositeACount && totalNumGPos >= minOppositeGCount){
								cytosinePassFlag = true;
								
							}
								
						}
						else if(totalNumPos >= 10 && ((double)(totalNumC + totalNumT)/(double)totalNumNeg) >= minCTContent){
							double freq = (double)(totalNumAPos)/(double)totalNumPos;
							if(freq <= minOppositeAFreq  && totalNumGPos >= minOppositeGCount){
								cytosinePassFlag = true;
							
							}
								
						}
						
					}
					else{
						if(totalNumNeg <= 0 && ((double)(totalNumC + totalNumT)/(double)totalNumPos) >= minCTContent){
							cytosinePassFlag = true;
						
							
						}
						else if(totalNumNeg > 0 && totalNumNeg < 10){
							if(totalNumANeg < minOppositeACount && totalNumGNeg >= minOppositeGCount){
								cytosinePassFlag = true;
								
							}
								
						}
						else if(totalNumNeg >= 10){
							double freq = (double)(totalNumANeg)/(double)totalNumNeg;
							if(freq <= minOppositeAFreq && totalNumGNeg >= minOppositeGCount){
								cytosinePassFlag = true;
							
							}
								
						}
					}
				}
				*/
				if(snpHomo && (tmpArray[5].equalsIgnoreCase("C") || tmpArray[5].equalsIgnoreCase("G"))){
					if(cytosineFlag){
						if(cytosinePassFlag){
							truePC++;
					//		cytosineWriter.println(line + "\tHOM\tHOM");
							
						}
						else{
							falseNC++;
					//		cytosineWriter.println(line + "\tHOM\tHET");
							System.err.println("FN: " + line);
						}	
					
					}
					else{
						if(cytosinePassFlag){
							truePC++;
							
					//		cytosineWriter.println(line + "\tHOM\tHOM");
							/*
							if(cpgPassFlag){
								truePCpg++;
								cpgWriter.println(line + "\tHOM\tHOM");
							}
							if(cphPassFlag){
								truePCph++;
								cphWriter.println(line + "\tHOM\tHOM");
							}*/
						}
						else{
							falseNC++;
							System.err.println("FN: " + line);
					//		cytosineWriter.println(line + "\tHOM\tHET");
							/*
							if(cpgPassFlag){
								falseNCpg++;
								cpgWriter.println(line + "\tHOM\tHET");
							}
							if(cphPassFlag){
								falseNCph++;
								cphWriter.println(line + "\tHOM\tHET");
							}*/
						}
						//homoCarrayNorefC++;
						//System.out.println("homo C in array, but no ref C: " + line + "\trefNegRecord: " + refNegRecord + "\trefPosRecord: " + refPosRecord + "\tCytosine negative: " + cytosineNegFlag);
					}
				}
				else if(!snpHomo && (tmpArray[5].equalsIgnoreCase("C") || tmpArray[5].equalsIgnoreCase("G") || tmpArray[6].equalsIgnoreCase("C") || tmpArray[6].equalsIgnoreCase("G"))){
					//if(cytosineFlag){
						if(cytosinePassFlag){
							falsePC++;
							System.err.println("FP: " + line);
					//		cytosineWriter.println(line + "\tHET\tHOM");
						/*
							if(cpgPassFlag){
								falsePCpg++;
								cpgWriter.println(line + "\tHET\tHOM");
							}
							if(cphPassFlag){
								falsePCph++;
								cphWriter.println(line + "\tHET\tHOM");
							}*/
							
						}
						else{
							trueNC++;
					//		cytosineWriter.println(line + "\tHET\tHET");
							/*
							if(cpgPassFlag){
								trueNCpg++;
								cpgWriter.println(line + "\tHET\tHET");
							}
							if(cphPassFlag){
								trueNCph++;
								cphWriter.println(line + "\tHET\tHET");
							}*/
						}	
					//}
					
				}
				
					
			}
			inputSam.close();
		//	cytosineWriter.close();
		
			double falsePositiveC = (double)falsePC/(double)snpHeterLociCafterCTreadcovFilter;
			double falseNegativeC = (double)falseNC/(double)snpHomoLociCafterCTreadcovFilter;
			int falsePCtrue = snpHeterLociC - trueNC;
			int falseNCtrue = snpHomoLociC - truePC;
			double falsePositiveCtrue = (double)falsePCtrue/(double)snpHeterLociC;
			double falseNegativeCtrue = (double)falseNCtrue/(double)snpHomoLociC;

			System.out.printf("Cytosine true positive number is: %d\n",truePC);
			System.out.printf("Cytosine true negative number is: %d\n",trueNC);
			System.out.printf("Cytosine false positive rate is: %f\n",falsePositiveC);
			System.out.printf("Cytosine false negative rate is: %f\n",falseNegativeC);
			System.out.printf("Cytosine false positive number is: %d\n",falsePC);
			System.out.printf("Cytosine false negative number is: %d\n",falseNC);
			System.out.printf("Cytosine snpHeterLociCafterCTreadcovFilter number is: %d\n",snpHeterLociCafterCTreadcovFilter);
			System.out.printf("Cytosine snpHomoLociCafterCTreadcovFilter number is: %d\n",snpHomoLociCafterCTreadcovFilter);
			
			System.out.printf("Cytosine true false positive rate is: %f\n",falsePositiveCtrue);
			System.out.printf("Cytosine true false negative rate is: %f\n",falseNegativeCtrue);
			System.out.printf("Cytosine false positive number is: %d\n",falsePCtrue);
			System.out.printf("Cytosine false negative number is: %d\n",falseNCtrue);
			System.out.printf("Cytosine snpHeterLociC number is: %d\n",snpHeterLociC);
			System.out.printf("Cytosine snpHomoLociC number is: %d\n",snpHomoLociC);
			
		//	System.out.printf("homo C in array, but no ref C: %d\n",homoCarrayNorefC);
			
			//System.out.printf("over");
	}

}
