package edu.usc.epigenome.scripts;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
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

public class SamToSnpAllNew {

	/**
	 * @param args
	 */
	final private static String prefix = "methylCGsRich_ASM_";
	final private static String USAGE = "SamToSnpAll [opts] sampleName file1.bam file2.bam ...";

	final private static int PURGE_INTERVAL = 20000; // We purge our stored Cpgs once we get this many bases past them.
	
	
	/**
	 * object vars
	 */

	
	/**
	 * @param args
	 */
	@Option(name="-chrom",multiValued=true,usage="One or more chroms, eg. --chrom chr1 --chrom chr5")
	protected List<String> chrs = new ArrayList<String>(25);
	@Option(name="-minMapQ",usage="minimum mapping quality (default 30)")
	protected int minMapQ = 10;
	@Option(name="-minBaseQual",usage="minimum Base quality (default 10)")
	protected int minBaseQual = 10;
	@Option(name="-minReadCov",usage="minimum read coverage (default 4)")
	protected int minReadCov = 4;
	@Option(name="-minAlleleCount",usage="minimum Allele Count (default 2)")
	protected static int minAlleleCount = 2;
	@Option(name="-minAlleleFreq",usage="minimum BAllele Frequency (default 0.2)")
	protected static double minAlleleFreq = 0.15;
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
		new SamToSnpAllNew().doMain(args);
	}

	public void doMain(String[] args)
	throws Exception {

		CmdLineParser parser = new CmdLineParser(this);
		parser.setUsageWidth(80);
		try
		{
			parser.parseArgument(args);
			if (stringArgs.size() < 2) throw new CmdLineException(USAGE);
			
			
		}
		catch (CmdLineException e)
		{
			System.err.println(e.getMessage());
			parser.printUsage(System.err);
			System.err.println();
			return;
		}

		String sampleName = stringArgs.remove(0);
		if (chrs.size()==0) chrs = MethylDbUtils.CHROMS;
		for (final String chr : chrs)
		{
			String tableName = prefix + sampleName + "_" + chr;

			for (final String fn : stringArgs)
			{
				File inputSamOrBamFile = new File(fn);

				SAMFileReader inputSam = new SAMFileReader(inputSamOrBamFile);
				inputSam.setValidationStringency(SAMFileReader.ValidationStringency.SILENT);
				boolean canPurge = ((inputSam.hasIndex()) && (stringArgs.size() == 1));
				Logger.getLogger(Logger.GLOBAL_LOGGER_NAME).info(String.format("Able to purge at interval %d? %s\n", 
						PURGE_INTERVAL, (canPurge) ? "Yes" : "Purging not available - either unsorted BAM or multiple BAMs for the same chrom"));
				
				TreeSet<Integer> allelePosition = allelePos( inputSam, canPurge, chr );
				inputSam.close();
				inputSam = new SAMFileReader(inputSamOrBamFile);
				inputSam.setValidationStringency(SAMFileReader.ValidationStringency.SILENT);
				outputSNP(allelePosition, tableName, chr, inputSam);
				inputSam.close();
				System.err.println("finished!");
			}
		}
	}

		protected TreeSet<Integer> allelePos( SAMFileReader inputSam, boolean canPurge, String chr)
		throws Exception{
			CloseableIterator<SAMRecord> chrItForAllele = inputSam.query(chr, 0, 0, false);
			int lastBaseSeen = 0;
			int lastPurge = 0;
			int recCounter = 0;
			TreeSet<Integer> allelePosition = new TreeSet<Integer>();
			
			
			recordSub: while (chrItForAllele.hasNext())
			{
				if (canPurge && (lastBaseSeen > (lastPurge+PURGE_INTERVAL)))
				{
					
					lastPurge = lastBaseSeen;
				}
				
				SAMRecord samRecord = chrItForAllele.next();

				// Filter low qual
				int mapQual = samRecord.getMappingQuality();
				byte[] baseQual = samRecord.getBaseQualities();
				boolean unmapped = samRecord.getReadUnmappedFlag();
				if (unmapped || (mapQual <= minMapQ))
				{
					continue recordSub;
				}


				String seq = PicardUtils.getReadString(samRecord, true);
				
				recCounter++;
				if ((recCounter % 1E5)==0)
				{
					System.err.printf("On new record #%d\n",recCounter);
					if (canPurge) System.gc();
				}


				try
				{
					String ref = PicardUtils.refStr(samRecord, true);

					if (seq.length() != ref.length())
					{
						System.err.println("SeqLen(" + seq.length() + ") != RefLen(" + ref.length() + ")");
						System.err.println(seq + "\n" + ref);
					}

					boolean negStrand = samRecord.getReadNegativeStrandFlag();
					int alignmentS = samRecord.getAlignmentStart();
					int	onRefCoord = (negStrand) ? samRecord.getUnclippedEnd() : alignmentS; 
					
					if (alignmentS < lastBaseSeen)
					{
						System.err.printf("BAM must be ordered in order: %d<%d\n",alignmentS, lastBaseSeen);
						System.exit(1);
					}
					lastBaseSeen = alignmentS;

					int seqLen = Math.min(seq.length(), ref.length());
					
					for (int i = 0; i < seqLen; i++){
						char refi = ref.charAt(i);
						byte baseQS = (negStrand) ? baseQual[seqLen-1-i] : baseQual[i];	
						if(baseQS > minBaseQual){
								if ((PicardUtils.isGuanine(i,ref) && PicardUtils.isAdenine(i,seq)) || (PicardUtils.isAdenine(i,ref) && PicardUtils.isGuanine(i,seq)) || (PicardUtils.isGuanine(i,ref) && PicardUtils.isCytosine(i, seq, true)) || (PicardUtils.isAdenine(i,ref) && PicardUtils.isCytosine(i, seq, true)) || (PicardUtils.isGuanine(i,seq) && PicardUtils.isCytosine(i, ref, true)) || (PicardUtils.isAdenine(i,seq) && PicardUtils.isCytosine(i, ref, true))){
									allelePosition.add(onRefCoord);
								}
						}
						if (refi == '-')
						{
							// It's a deletion in reference, don't advance
							allelePosition.remove(onRefCoord);
							allelePosition.remove(onRefCoord-1);
							allelePosition.remove(onRefCoord-2);
							allelePosition.remove(onRefCoord-3);
							int inc = (negStrand) ? -3 : 3;
							onRefCoord += inc;
							i+=3;
							
						}
						else
						{
							int inc = (negStrand) ? -1 : 1;
							onRefCoord += inc;
						}
					}
				}
				catch (Exception f)
				{
					System.err.println("-----------------------------------------");
					System.err.println("Couldn't handle seq #" + recCounter);
					System.err.println(seq);
					f.printStackTrace(System.err);
					System.err.println("-----------------------------------------");
				}
			}
			chrItForAllele.close();
			return allelePosition;
			
		}
	
		
		protected void outputSNP(TreeSet<Integer> allelePosition, String tableName, String chr, SAMFileReader inputSam)
		throws FileNotFoundException{
			Iterator<Integer> it = allelePosition.iterator();
			System.err.println("-----------------------------------------");
			String fn = tableName + "_SNP_all_afterBaseQfilter" + ".txt";
			PrintWriter writer = new PrintWriter(new File(fn));
			int count = 0;
			//boolean xFlag = false;
			//if(chr.equalsIgnoreCase("chrX")){
			//	xFlag = true;
			//}
			
			while(it.hasNext()){	
				Integer snpPosition = it.next();
				boolean agFlag = false;
				boolean agNegFlag = false;
				
				int totalNum = 0;
				int alleleNum = 0;
				int totalNumPos = 0;
				int alleleNumPos = 0;
				int totalNumNeg = 0;
				int alleleNumNeg = 0;
				//char refBase;
				//char alleleBase;
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
						if((PicardUtils.isGuanine(i,ref) && PicardUtils.isAdenine(i,seq)) || (PicardUtils.isAdenine(i,ref) && PicardUtils.isGuanine(i,seq))){
							agFlag = true;
							//refBase = ref.charAt(i);
							//alleleBase = seq.charAt(i);
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
							//refBase = ref.charAt(i);
							//alleleBase = seq.charAt(i);
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
				if(totalNum < minReadCov){
						continue;
					}
				//if(totalNum < minReadCov & (!xFlag || female)){
				//	continue;
				//}
				//else if(totalNum < minReadCov/2 & (xFlag && !female)){
				//	continue;
				//}
				
				double freq1 = (double)alleleNum/(double)totalNum;
				double freq2 = (double)(totalNum-alleleNum)/(double)totalNum;
				if(totalNum >10 && freq1 >= minAlleleFreq && freq2 >= minAlleleFreq && alleleNum >= minAlleleCount && (totalNum - alleleNum) >= minAlleleCount){
					writer.printf("%d\t%d\t%d\n",snpPosition,totalNum,alleleNum);
					//System.out.println(snpPosition);
					count++;
				 }
				else if(totalNum <= 10 && alleleNum >= minAlleleCount && (totalNum - alleleNum) >= minAlleleCount){
					writer.printf("%d\t%d\t%d\n",snpPosition,totalNum,alleleNum);
					//System.out.println(snpPosition);
					count++;
				}
				
			}
			System.out.println(count);
			writer.close();
		}
		
		
		
}
