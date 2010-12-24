package edu.usc.epigenome.scripts;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.TreeMap;
import java.util.TreeSet;

import java.util.logging.Logger;

import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;

import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.util.CloseableIterator;
import edu.usc.epigenome.genomeLibs.PicardUtils;
import edu.usc.epigenome.genomeLibs.MethylDb.MethylDbUtils;

public class SamToFpFn {

	/**
	 * @param args
	 */
		
		//final private static String prefix = "methylCGsRich_ASM_AllSnp_";
		final private static String USAGE = "SamToBaseQ [opts] sampleName snpFileName bamFilePath ...";

		final private static int PURGE_INTERVAL = 20000; // We purge our stored Cpgs once we get this many bases past them.
		//final private static int ALLELE_GA_NUMBER = 1; //for each reads, when there are more than 1 GA position is different from reference sequence, we define it belongs to another allele. 
		
		/**
		 * @param args
		 */
		
		@Option(name="-minMapQ",usage="minimum mapping quality (default 30)")
		protected int minMapQ = 30;
		@Option(name="-minReadCov",usage="minimum read coverage (default 4)")
		protected int minReadCov = 4;
		@Option(name="-minBaseQual",usage="minimum Base quality (default 10)")
		protected int minBaseQual = 20;
		@Option(name="-minAlleleCount",usage="minimum Allele Count (default 2)")
		protected static int minAlleleCount = 1;
		@Option(name="-minAlleleFreq",usage="minimum BAllele Frequency (default 0.2)")
		protected static double minAlleleFreq = 0.10;
		@Option(name="-female",usage="the sample is female or male (default false for male)")
		protected boolean female = false;
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
			new SamToFpFn().doMain(args);
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
			String snpFileName = stringArgs.remove(0);
			String bamFilePath = stringArgs.remove(0);
			BufferedReader br = new BufferedReader(new FileReader(snpFileName));
			String line;
			int heterLoci = 0;
			int homoLoci = 0;
			int snpHeterLoci = 0;
			int snpHomoLoci = 0;
			int trueP = 0;
			int falseP = 0;
			int trueN = 0;
			int falseN = 0;
			int truePX = 0;
			int falsePX = 0;
			int trueNX = 0;
			int falseNX = 0;
			int heterLociX = 0;
			int homoLociX = 0;
			int snpHeterLociX = 0;
			int snpHomoLociX = 0;
			int snpHeterLociAll = 0;
			String preChr = "chr1";
			String fn = sampleName + "Merge_chr1.NODUPS.sorted.calmd.NODUPS.bam";
			File inputSamOrBamFile = new File(bamFilePath,fn);
			SAMFileReader inputSam = new SAMFileReader(inputSamOrBamFile);
			inputSam.setValidationStringency(SAMFileReader.ValidationStringency.SILENT);
			String fnAllSnp = sampleName + "_SNP.filterCNV.filterSeqDepth." + minMapQ + "-" + minBaseQual + "-" + minReadCov + "-" + minAlleleCount + "-" + minAlleleFreq + "." + "txt";
			String fnFalseSnp = sampleName + "_SNP.filterCNV.filterSeqDepth.falseSnp." + minMapQ + "-" + minBaseQual + "-" + minReadCov + "-" + minAlleleCount + "-" + minAlleleFreq + "." + "txt";
			PrintWriter allSnpWriter = new PrintWriter(new File(fnAllSnp));
			PrintWriter falseSnpWriter = new PrintWriter(new File(fnFalseSnp));
			
			while( (line = br.readLine()) != null){
				String[] tmpArray = line.split("\t");
				String chr = "chr" + tmpArray[1];
				boolean xFlag = false;
				if(chr.equalsIgnoreCase("chrX")){
					xFlag = true;
				}
				if(!chr.equalsIgnoreCase(preChr)){
					inputSam.close();
					fn = sampleName + "Merge_" + chr + ".NODUPS.sorted.calmd.NODUPS.bam";
					File inputNewSamOrBamFile = new File(bamFilePath,fn);
					System.err.println(inputNewSamOrBamFile.getName());
					inputSam = new SAMFileReader(inputNewSamOrBamFile);
					inputSam.setValidationStringency(SAMFileReader.ValidationStringency.SILENT);
					preChr = chr;
				}
				Integer snpPosition =(int) Math.round(Double.parseDouble(tmpArray[2]));
				//System.err.println(Double.parseDouble(tmpArray[2]));
				//System.err.println(snpPosition);
				boolean snpHomo = tmpArray[3].equalsIgnoreCase(tmpArray[4]);
				boolean agFlag = false;
				boolean agNegFlag = false;
				
				int totalNum = 0;
				int alleleNum = 0;
				int totalNumPos = 0;
				int alleleNumPos = 0;
				int totalNumNeg = 0;
				int alleleNumNeg = 0;
				String posSeq = "-";
				String negSeq = "-";
				
				if(!snpHomo){
					snpHeterLociAll++;
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
						if(i >= 1 && i < seqLen-1){
							if(negStrand){
								negSeq = ref.substring(i-1,i+2);
							}
							else{
								posSeq = ref.substring(i-1,i+2);
							}
						}
						totalNum++;
						if(negStrand){
							totalNumNeg++;
						}
						else{
							totalNumPos++;
						}
						if((PicardUtils.isGuanine(i,ref) && PicardUtils.isAdenine(i,seq)) || (PicardUtils.isAdenine(i,ref) && PicardUtils.isGuanine(i,seq))){
							agFlag = true;
							if(negStrand){
								alleleNumNeg++;
								agNegFlag = true;
							}
							else{
								alleleNumPos++;
							}
						}
						else if ((PicardUtils.isGuanine(i,ref) && PicardUtils.isCytosine(i, seq, true)) || (PicardUtils.isAdenine(i,ref) && PicardUtils.isCytosine(i, seq, true)) || (PicardUtils.isGuanine(i,seq) && PicardUtils.isCytosine(i, ref, true)) || (PicardUtils.isAdenine(i,seq) && PicardUtils.isCytosine(i, ref, true))){	
							//alleleNum++;// allele reads number
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
				if(totalNum < minReadCov & (!xFlag || female)){
					continue;
				}
				else if(totalNum < minReadCov/2 & (xFlag && !female)){
					continue;
				}
					/*if(totalNumPos > 0){
						double freq1Pos = (double)alleleNumPos/(double)totalNumPos;
						double freq2Pos = (double)(totalNumPos-alleleNumPos)/(double)totalNumPos;
						if(totalNumPos >10 && freq1Pos >= minAlleleFreq && freq2Pos >= minAlleleFreq && alleleNumPos >= minAlleleCount && (totalNumPos - alleleNumPos) >= minAlleleCount){
							totalNum = totalNumPos;
							alleleNum = alleleNumPos;
						}
						else if(totalNumPos <= 10 && alleleNumPos >= minAlleleCount && (totalNumPos - alleleNumPos) >= minAlleleCount){
							totalNum = totalNumPos;
							alleleNum = alleleNumPos;
						}
						else{
							
						}
					}
					if(totalNumNeg > 0){
						double freq1Neg = (double)alleleNumNeg/(double)totalNumNeg;
						double freq2Neg = (double)(totalNumNeg-alleleNumNeg)/(double)totalNumNeg;
						if(totalNumNeg >10 && freq1Neg >= minAlleleFreq && freq2Neg >= minAlleleFreq && alleleNumNeg >= minAlleleCount && (totalNumNeg - alleleNumNeg) >= minAlleleCount){
							totalNum = totalNumNeg;
							alleleNum = alleleNumNeg;
						}
						else if(totalNumNeg <= 10 && alleleNumNeg >= minAlleleCount && (totalNumNeg - alleleNumNeg) >= minAlleleCount){
							totalNum = totalNumNeg;
							alleleNum = alleleNumNeg;
						}
						else{
							
						}
					}*/
					
					double freq1 = (double)alleleNum/(double)totalNum;
					double freq2 = (double)(totalNum-alleleNum)/(double)totalNum;
					
					allSnpWriter.printf("%s\t%s\t%s\t%d\t%d\t%d\t%d\n",line,posSeq,negSeq,totalNumPos,alleleNumPos,totalNumNeg,alleleNumNeg);
					if(snpHomo){
						snpHomoLoci++;
						if(xFlag){
							snpHomoLociX++;
						}
					}
					else{
						snpHeterLoci++;
						if(xFlag){
							snpHeterLociX++;
						}
					}
					boolean seqHomo;
					
					if(totalNum >10 && freq1 >= minAlleleFreq && freq2 >= minAlleleFreq && alleleNum >= minAlleleCount && (totalNum - alleleNum) >= minAlleleCount){
						seqHomo = false;
						heterLoci++;
						if(xFlag){
							heterLociX++;
						}
					 }
					else if(totalNum <= 10 && alleleNum >= minAlleleCount && (totalNum - alleleNum) >= minAlleleCount){
						seqHomo = false;
						heterLoci++;
						if(xFlag){
							heterLociX++;
						}
					}
					else{
						seqHomo = true;
						homoLoci++;
						if(xFlag){
							homoLociX++;
						}
					}
					if(snpHomo && seqHomo){
						trueN++;
						if(xFlag){
							trueNX++;
						}
					}
					else if(!snpHomo && !seqHomo){
						trueP++;
						if(xFlag){
							truePX++;
						}
					}
					else if(snpHomo && !seqHomo){
						falseP++;
						if(xFlag){
							falsePX++;
						}

							System.err.println(line);
							falseSnpWriter.printf("%s\t%s\t%s\t%d\t%d\t%d\t%d\n",line,posSeq,negSeq,totalNumPos,alleleNumPos,totalNumNeg,alleleNumNeg);
					}
					else if(!snpHomo && seqHomo){
						falseN++;
						if(xFlag){
							falseNX++;
						}
						System.err.println(line);
						falseSnpWriter.printf("%s\t%s\t%s\t%d\t%d\t%d\t%d\n",line,posSeq,negSeq,totalNumPos,alleleNumPos,totalNumNeg,alleleNumNeg);
						if(debug){
							System.err.println(line);
							
							System.err.println(alleleNum);
							System.err.println(totalNum);
							System.err.println(!seqHomo);
							System.err.println(!snpHomo);
							System.err.println(heterLoci);
							System.err.println(freq1);
							System.err.println(freq2);
							System.err.println(alleleNumPos);
							System.err.println(totalNumPos);
							System.err.println(alleleNumNeg);
							System.err.println(totalNumNeg);
						}
					}
					/*if(!snpHomo){
						System.err.println(line);
						System.err.println(alleleNum);
						System.err.println(totalNum);
						System.err.println(!seqHomo);
						System.err.println(!snpHomo);
						System.err.println(heterLoci);
						System.err.println(freq1);
						System.err.println(freq2);
					}*/
					
			}
			inputSam.close();
			allSnpWriter.close();
			falseSnpWriter.close();
			double falsePositive = (double)falseP/(double)heterLoci;
			double falseNegative = (double)falseN/(double)snpHeterLoci;
			double FN = (double)falseN/(double)snpHeterLociAll;
			System.out.printf("heterozygous loci in 1M SNP array is: %d\n",snpHeterLoci);
			System.out.printf("homozygous loci in 1M SNP array is: %d\n",snpHomoLoci);
			System.out.printf("heterozygous loci in sequencing is: %d\n",heterLoci);
			System.out.printf("homozygous loci in sequencing is: %d\n",homoLoci);
			System.out.printf("true positive number is: %d\n",trueP);
			System.out.printf("true negative number is: %d\n",trueN);
			System.out.printf("false positive rate is: %f\n",falsePositive);
			System.out.printf("false negative rate is: %f\n",falseNegative);
			System.out.printf("FN rate is: %f\n",FN);
			System.out.printf("X chromosome heterozygous locis in 1M SNP array is: %d\n",snpHeterLociX);
			System.out.printf("X chromosome homozygous locis in 1M SNP array is: %d\n",snpHomoLociX);
			System.out.printf("X chromosome heterozygous locis in sequencing is: %d\n",heterLociX);
			System.out.printf("X chromosome homozygous locis in sequencing is: %d\n",homoLociX);
			System.out.printf("X chromosome true positive numbers is: %d\n",truePX);
			System.out.printf("X chromosome true negative numbers is: %d\n",trueNX);
	}
}
