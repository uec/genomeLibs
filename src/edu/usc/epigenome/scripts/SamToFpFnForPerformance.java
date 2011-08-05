package edu.usc.epigenome.scripts;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
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

public class SamToFpFnForPerformance {
	
		
		//final private static String prefix = "methylCGsRich_ASM_AllSnp_";
		final private static String USAGE = "SamToFpFn [opts] snpFileName bamFile ...";

		final private static int PURGE_INTERVAL = 20000; // We purge our stored Cpgs once we get this many bases past them.
		//final private static int ALLELE_GA_NUMBER = 1; //for each reads, when there are more than 1 GA position is different from reference sequence, we define it belongs to another allele. 
		
		/**
		 * @param args
		 */
		
		@Option(name="-minMapQ",usage="minimum mapping quality (default 30)")
		protected int minMapQ = 30;
		@Option(name="-minReadCov",usage="minimum read coverage (default 1)")
		protected int minReadCov = 1;
		@Option(name="-minBaseQual",usage="minimum Base quality (default 0)")
		protected int minBaseQual = 0;
		@Option(name="-minAlleleCount",usage="minimum Allele Count (default 2)")
		protected static double minAlleleCount = 2;
		@Option(name="-minAlleleFreq",usage="minimum BAllele Frequency (default 0.2)")
		protected static double minAlleleFreq = 0.2;
		//@Option(name="-female",usage="the sample is female or male (default false for male)")
		//protected boolean female = false;
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
			SamToFpFnForPerformance test = new SamToFpFnForPerformance();
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
			int seqHeterLoci = 0;
			int seqHomoLoci = 0;
			int snpHeterLoci = 0;
			int snpHomoLoci = 0;
			int trueP = 0;
			int falseP = 0;
			int trueN = 0;
			int falseN = 0;
			
			int recCounter = 1;
			//int snpHeterLociAll = 0;
			//String fn = sampleName + "Merge_chr1.NODUPS.sorted.calmd.NODUPS.bam";
			File inputSamOrBamFile = new File(bamFile);
			SAMFileReader inputSam = new SAMFileReader(inputSamOrBamFile);
			inputSam.setValidationStringency(SAMFileReader.ValidationStringency.SILENT);
			//String fnAllSnp = sampleName + "_SNP.filterCNV.filterSeqDepth." + minMapQ + "-" + minBaseQual + "-" + minReadCov + "-" + minAlleleCount + "-" + minAlleleFreq + "." + "txt";
			//String fnFalseSnp = sampleName + "_SNP.filterCNV.filterSeqDepth.falseSnp." + minMapQ + "-" + minBaseQual + "-" + minReadCov + "-" + minAlleleCount + "-" + minAlleleFreq + "." + "txt";
			//PrintWriter allSnpWriter = new PrintWriter(new File(fnAllSnp));
			//PrintWriter falseSnpWriter = new PrintWriter(new File(fnFalseSnp));
			
			while( (line = br.readLine()) != null){
				String[] tmpArray = line.split("\t");
				String chr = "chr" + tmpArray[1];
				
				if ((recCounter % 10000)==0)
				{
					System.err.printf("On new record #%d\n",recCounter);
					inputSam.close();
					
					System.gc();
					inputSam = new SAMFileReader(inputSamOrBamFile);
					inputSam.setValidationStringency(SAMFileReader.ValidationStringency.SILENT);
					
				}
				
				Integer snpPosition =(int) Math.round(Double.parseDouble(tmpArray[2]));
				//System.err.println(Double.parseDouble(tmpArray[2]));
				//System.err.println(snpPosition);
				boolean snpHomo = tmpArray[3].equalsIgnoreCase(tmpArray[4]);
				boolean seqHomo = true;
					
				boolean agFlag = false;
				boolean agNegFlag = false;
				boolean noAgFlag = false;
				// this is just for cheat, to see the sepcific sequence depth filter to AG/TC SNP would make this kind of SNP to normal number...(normal number is 4X than other SNP) 
				//if((tmpArray[5].equalsIgnoreCase("A") && tmpArray[6].equalsIgnoreCase("G")) || (tmpArray[5].equalsIgnoreCase("T") && tmpArray[6].equalsIgnoreCase("C"))){
				//	agFlag = true;
					//System.err.println(line);
				//}
				int totalNum = 0;
				int alleleNum = 0;
				int totalNumPos = 0;
				int alleleNumPos = 0;
				int totalNumNeg = 0;
				int alleleNumNeg = 0;
				
				
				if(snpHomo){
					snpHomoLoci++;
					
				}
				else{
					snpHeterLoci++;
					
				}
				
				CloseableIterator<SAMRecord> chrIt = inputSam.queryOverlapping(chr, snpPosition, snpPosition);
				record: while (chrIt.hasNext())
				{	
					SAMRecord samRecord = chrIt.next();
					int mapQual = samRecord.getMappingQuality();
					byte[] baseQual = samRecord.getBaseQualities();
					boolean unmapped = samRecord.getReadUnmappedFlag();
					if (unmapped || (mapQual <= minMapQ))
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
						if( baseQS <= minBaseQual )
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
						if(PicardUtils.isCytosine(i, ref, true) && PicardUtils.isCytosine(i, seq, true) && !noAgFlag){
							agFlag = true;
							if(!negStrand){
								agNegFlag = true;
							}
						}
						if((PicardUtils.isGuanine(i,ref) && seq.charAt(i) == 'A') || (ref.charAt(i) == 'A' && PicardUtils.isGuanine(i,seq))){
							agFlag = true;
							noAgFlag = false;
							
							if(negStrand){
								alleleNumNeg++;
								agNegFlag = true;
							}
							else{
								alleleNumPos++;
							}
						}
						else if ((PicardUtils.isGuanine(i,ref) && PicardUtils.isCytosine(i, seq, true)) || (ref.charAt(i) == 'A' && PicardUtils.isCytosine(i, seq, true)) || (PicardUtils.isGuanine(i,seq) && PicardUtils.isCytosine(i, ref, true)) || (seq.charAt(i) == 'A' && PicardUtils.isCytosine(i, ref, true))){	
							
							noAgFlag = true;
							if(negStrand){
								alleleNumNeg++;
							}
							else{
								alleleNumPos++;
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
				if(agFlag){
					if(agNegFlag){
						alleleNum = alleleNumNeg;
						totalNum = totalNumNeg;
					}
					else{
						alleleNum = alleleNumPos;
						totalNum = totalNumPos;
					}
				}
				else{
					alleleNum = alleleNumPos + alleleNumNeg;
					totalNum = totalNumPos + totalNumNeg;
				}
				if(totalNum >= minReadCov ){
					double freq1 = (double)alleleNum/(double)totalNum;
					double freq2 = (double)(totalNum-alleleNum)/(double)totalNum;
					if(totalNum >10 && freq1 >= minAlleleFreq && freq2 >= minAlleleFreq && alleleNum >= minAlleleCount && (totalNum - alleleNum) >= minAlleleCount){
						seqHomo = false;
					 }
					else if(totalNum <= 10 && alleleNum >= minAlleleCount && (totalNum - alleleNum) >= minAlleleCount){
						seqHomo = false;
					}
					else{
						seqHomo = true;
					}
				}
				
				if(seqHomo){
					seqHomoLoci++;
					
				}
				else{
					seqHeterLoci++;
					
				}
				
					if(snpHomo && seqHomo){
						trueN++;
						
					}
					else if(!snpHomo && !seqHomo){
						trueP++;
						
					}
					else if(snpHomo && !seqHomo){
						falseP++;
						
						System.err.println(line);
							
					}
					else if(!snpHomo && seqHomo){
						falseN++;
						
						System.err.println(line);
						
						if(debug){
							System.err.println(line);
							
							System.err.println(alleleNum);
							System.err.println(totalNum);
							System.err.println(!seqHomo);
							System.err.println(!snpHomo);
							System.err.println(seqHeterLoci);
							//System.err.println(freq1);
							//System.err.println(freq2);
							System.err.println(alleleNumPos);
							System.err.println(totalNumPos);
							System.err.println(alleleNumNeg);
							System.err.println(totalNumNeg);
						}
					}
					
					
			}
			inputSam.close();
			//allSnpWriter.close();
			//falseSnpWriter.close();
			double falsePositive = (double)falseP/(double)snpHomoLoci;
			double falseNegative = (double)falseN/(double)snpHeterLoci;

			System.out.printf("heterozygous loci in 1M SNP array is: %d\n",snpHeterLoci);
			System.out.printf("homozygous loci in 1M SNP array is: %d\n",snpHomoLoci);
			System.out.printf("heterozygous loci in sequencing is: %d\n",seqHeterLoci);
			System.out.printf("homozygous loci in sequencing is: %d\n",seqHomoLoci);
			System.out.printf("true positive number is: %d\n",trueP);
			System.out.printf("true negative number is: %d\n",trueN);
			System.out.printf("false positive rate is: %f\n",falsePositive);
			System.out.printf("false negative rate is: %f\n",falseNegative);
			
			//System.out.printf("over");
	}

}
