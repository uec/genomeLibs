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

public class SamToFpFnMethyByBismark {

	/**
	 * @param args
	 */
		
		//final private static String prefix = "methylCGsRich_ASM_AllSnp_";
		final private static String USAGE = "SamToFpFnMethyByBismark [opts] snpFileName bamFile ...";

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
		@Option(name="-minCreadsCount",usage="minimum cytosine Count (default 0)")
		protected static double minCreadsCount = 0;
		
		//@Option(name="-cpg",usage=" seperately consider CpG and CpH sites (default false)")
		//protected boolean cpg = false;
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
			SamToFpFnMethyByBismark test = new SamToFpFnMethyByBismark();
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
			
			
			//int snpHeterLociAll = 0;
			int recCounter = 1;
			//String fn = sampleName + "Merge_chr1.NODUPS.sorted.calmd.NODUPS.bam";
			File inputSamOrBamFile = new File(bamFile);
			SAMFileReader inputSam = new SAMFileReader(inputSamOrBamFile);
			inputSam.setValidationStringency(SAMFileReader.ValidationStringency.SILENT);
			//String fnAllSnp = sampleName + "_SNP.filterCNV.filterSeqDepth." + minMapQ + "-" + minBaseQual + "-" + minReadCov + "-" + minAlleleCount + "-" + minAlleleFreq + "." + "txt";
			//String fnFalseSnp = sampleName + "_SNP.filterCNV.filterSeqDepth.falseSnp." + minMapQ + "-" + minBaseQual + "-" + minReadCov + "-" + minAlleleCount + "-" + minAlleleFreq + "." + "txt";
			//PrintWriter allSnpWriter = new PrintWriter(new File(fnAllSnp));
			//PrintWriter falseSnpWriter = new PrintWriter(new File(fnFalseSnp));
			snpFileName = snpFileName + ".minCreadsCount-" + minCreadsCount + ".bismark.txt";
			PrintWriter snpWriter = new PrintWriter(new File(snpFileName));
			
			
			while( (line = br.readLine()) != null){
				String[] tmpArray = line.split("\t");
				if((tmpArray[5].equalsIgnoreCase("A") && tmpArray[6].equalsIgnoreCase("A")) || (tmpArray[5].equalsIgnoreCase("T") && tmpArray[6].equalsIgnoreCase("T")) ||  (tmpArray[5].equalsIgnoreCase("A") && tmpArray[6].equalsIgnoreCase("T")) || (tmpArray[5].equalsIgnoreCase("T") && tmpArray[6].equalsIgnoreCase("A")))
					continue;
					
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
				
				boolean snpHomo = tmpArray[3].equalsIgnoreCase(tmpArray[4]);
				boolean seqHomo = false;
				
				
				int totalNum = 0;
				int totalNumNeg = 0;
				int totalNumPos = 0;
				int CNumPos = 0;
				int TNumPos = 0;
				int CNumNeg = 0;
				int TNumNeg = 0;
				char refSeq = 'N';
				boolean cytosineFlag = false;
				
				
				
				
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
						refSeq = ref.charAt(i);
						totalNum++;
						if(negStrand){
							totalNumNeg++;
						}
						else{
							totalNumPos++;
							
						}

						if(PicardUtils.isCytosine(i, ref, false) ){
							cytosineFlag = true;
							if(PicardUtils.isCytosine(i, seq, false)){
								if(!negStrand){
									CNumPos++;
								}
								else{
									CNumNeg++;
								}
							}
							else if(PicardUtils.isCytosine(i, seq, true)){
								if(!negStrand){
									TNumPos++;
								}
								else{
									TNumNeg++;
								}
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
				
				
					if(snpHomo){
						snpHomoLoci++;
						
					}
					else{
						snpHeterLoci++;
						
					}
				
					
				if(totalNum >= minReadCov ){
					if((CNumPos + CNumNeg) >= minCreadsCount && cytosineFlag){
						seqHomo = true;
						seqHomoLoci++;
						snpWriter.println(chr + "\t" + tmpArray[2] + "\t" + CNumPos + "\t" + TNumPos + "\t" + totalNumPos + "\t" + CNumNeg + "\t" + TNumNeg + "\t" + totalNumNeg);
					}
					else{
						seqHomo = false;
						seqHeterLoci++;
					}
				}
				else{
					continue;
				}

				if(seqHomo){
					seqHomoLoci++;
						
				}
				else{
					seqHeterLoci++;
						
				}

					if(snpHomo && seqHomo){
						trueP++;
		
					}
					else if(!snpHomo && !seqHomo){
						trueN++;
					
					}
					else if(snpHomo && !seqHomo){
						falseN++;
						
						System.err.println(line + "\tref:" + refSeq);
						
					}
					else if(!snpHomo && seqHomo){
						falseP++;
						
						System.err.println(line + "\tref:" + refSeq);
						
					}
					
					
					
					
			}
			inputSam.close();
			snpWriter.close();
			//allSnpWriter.close();
			//falseSnpWriter.close();
			double falsePositive = (double)falseP/(double)snpHeterLoci;
			double falseNegative = (double)falseN/(double)snpHomoLoci;
			

			System.out.printf("heterozygous C loci in 1M SNP array is: %d\n",snpHeterLoci);
			System.out.printf("homozygous C loci in 1M SNP array is: %d\n",snpHomoLoci);
			System.out.printf("heterozygous C loci in sequencing is: %d\n",seqHeterLoci);
			System.out.printf("homozygous C loci in sequencing is: %d\n",seqHomoLoci);
			System.out.printf("true positive number is: %d\n",trueP);
			System.out.printf("true negative number is: %d\n",trueN);
			System.out.printf("false positive rate is: %f\n",falsePositive);
			System.out.printf("false negative rate is: %f\n",falseNegative);
			System.out.printf("false positive number is: %d\n",falseP);
			System.out.printf("false negative number is: %d\n",falseN);
			
			
			//System.out.printf("over");
	}

}
