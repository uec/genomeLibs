package edu.usc.epigenome.scripts;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.TreeMap;
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

public class SamToSnpAll {

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
	protected int minMapQ = 30;
	@Option(name="-minBaseQual",usage="minimum Base quality (default 10)")
	protected int minBaseQual = 10;
	@Option(name="-minAlleleCount",usage="minimum Allele Count (default 3)")
	protected static int minAlleleCount = 3;
	@Option(name="-minAlleleFreq",usage="minimum BAllele Frequency (default 0.2)")
	protected static double minAlleleFreq = 0.20;
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
		new SamToSnpAll().doMain(args);
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
		if (chrs.size()==0) chrs = MethylDbUtils.CHROMS;
		for (final String chr : chrs)
		{
			String tableName = prefix + sampleName + "_" + chr;

			for (final String fn : stringArgs)
			{
				File inputSamOrBamFile = new File(fn);

				final SAMFileReader inputSam = new SAMFileReader(inputSamOrBamFile);
				inputSam.setValidationStringency(SAMFileReader.ValidationStringency.SILENT);
				
				
				
				
				// We can only purge if we are the only input file and we are sorted.
				boolean canPurge = ((inputSam.hasIndex()) && (stringArgs.size() == 1));
				//boolean canPurge = false;
				Logger.getLogger(Logger.GLOBAL_LOGGER_NAME).info(String.format("Able to purge at interval %d? %s\n", 
						PURGE_INTERVAL, (canPurge) ? "Yes" : "Purging not available - either unsorted BAM or multiple BAMs for the same chrom"));
				
				TreeMap<Integer,Integer[]> allelePosition = allelePos( inputSam, canPurge, chr );
				//TreeMap<Integer,List<Byte>> alleleReadsMem = alleleReads( inputSam, allelePosition, canPurge, chr );
				//outputSNP(allelePosition, alleleReadsMem, tableName);
				outputSNP(allelePosition, tableName);
				System.err.println("finished!");
			}
		}
	}

		protected TreeMap<Integer,Integer[]> allelePos( SAMFileReader inputSam, boolean canPurge, String chr)
		throws Exception{
			CloseableIterator<SAMRecord> chrItForAllele = inputSam.query(chr, 0, 0, false);
			int lastBaseSeen = 0;
			int lastPurge = 0;
			int recCounter = 0;
			
			TreeMap<Integer, Integer> readsDepth = new TreeMap<Integer, Integer>();
			//int j = 0;
			TreeMap<Integer,Integer[]> allelePosition = new TreeMap<Integer,Integer[]>();
			
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
				if (unmapped || (mapQual < minMapQ))
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
					//System.err.println(seq + "\n" + ref);

					boolean negStrand = samRecord.getReadNegativeStrandFlag();
					int alignmentS = samRecord.getAlignmentStart();
					int	onRefCoord = (negStrand) ? samRecord.getUnclippedEnd() : alignmentS; 
					//int readsStart = samRecord.getUnclippedEnd() <= alignmentS ? samRecord.getUnclippedEnd() : alignmentS;
					//int readsEnd = samRecord.getUnclippedEnd() <= alignmentS ? alignmentS : samRecord.getUnclippedEnd();
					if ((recCounter % 1E5)==0 && onRefCoord-10000 > 0)
					{
						System.err.printf("On new onRefCoord #%d\n",onRefCoord-10000);
						readsDepth.headMap(onRefCoord-10000).clear();
						
					}
					if (alignmentS < lastBaseSeen)
					{
						System.err.printf("BAM must be ordered in order: %d<%d\n",alignmentS, lastBaseSeen);
						System.exit(1);
					}
					lastBaseSeen = alignmentS;

					int seqLen = Math.min(seq.length(), ref.length());
					
					for (int i = 0; i < seqLen; i++){
						char refi = ref.charAt(i);
						//char seqi = seq.charAt(i);
						//if( onRefCoord == 2718538 && (PicardUtils.isAdenine(i,seq) || PicardUtils.isGuanine(i,seq))){
						//	System.err.println(baseQual[i]);
						//	System.err.println(i);
						//}
						byte baseQS = (negStrand) ? baseQual[seqLen-1-i] : baseQual[i];	
						Integer strandRefCoord = (negStrand) ? 0-onRefCoord : onRefCoord;
						
						if(baseQS > minBaseQual){
							if(allelePosition.containsKey(strandRefCoord)){
								Integer[] tempInt = allelePosition.get(strandRefCoord);
								//if (PicardUtils.isAdenine(i,seq) || PicardUtils.isGuanine(i,seq)){
									tempInt[1]++;//  total reads number
								//}
									
								//if ((PicardUtils.isGuanine(i,ref) && PicardUtils.isCytosine(i, seq, true)) || (PicardUtils.isAdenine(i,ref) && PicardUtils.isCytosine(i, seq, true))){	
								if ((PicardUtils.isGuanine(i,ref) && PicardUtils.isAdenine(i,seq)) || (PicardUtils.isAdenine(i,ref) && PicardUtils.isGuanine(i,seq)) || (PicardUtils.isGuanine(i,ref) && PicardUtils.isCytosine(i, seq, true)) || (PicardUtils.isAdenine(i,ref) && PicardUtils.isCytosine(i, seq, true)) || (PicardUtils.isGuanine(i,seq) && PicardUtils.isCytosine(i, ref, true)) || (PicardUtils.isAdenine(i,seq) && PicardUtils.isCytosine(i, ref, true))){
								//if ((PicardUtils.isGuanine(i,ref) && PicardUtils.isAdenine(i,seq)) || (PicardUtils.isAdenine(i,ref) && PicardUtils.isGuanine(i,seq))){	
									tempInt[0]++;// allele reads number
								}
								allelePosition.put(strandRefCoord,tempInt);

									
							}
							else{
								Integer[] tempInt = new Integer[2];
								//if (PicardUtils.isAdenine(i,seq) || PicardUtils.isGuanine(i,seq)){
									tempInt[1] = readsDepth.get(strandRefCoord);
									if(tempInt[1] != null)
										tempInt[1]++;
									else
										tempInt[1] = 1;
									readsDepth.put(strandRefCoord,tempInt[1]);
								//}
								
									//if ((PicardUtils.isGuanine(i,ref) && PicardUtils.isCytosine(i, seq, true)) || (PicardUtils.isAdenine(i,ref) && PicardUtils.isCytosine(i, seq, true))){	
								if ((PicardUtils.isGuanine(i,ref) && PicardUtils.isAdenine(i,seq)) || (PicardUtils.isAdenine(i,ref) && PicardUtils.isGuanine(i,seq)) || (PicardUtils.isGuanine(i,ref) && PicardUtils.isCytosine(i, seq, true)) || (PicardUtils.isAdenine(i,ref) && PicardUtils.isCytosine(i, seq, true)) || (PicardUtils.isGuanine(i,seq) && PicardUtils.isCytosine(i, ref, true)) || (PicardUtils.isAdenine(i,seq) && PicardUtils.isCytosine(i, ref, true))){	
								//if ((PicardUtils.isGuanine(i,ref) && PicardUtils.isAdenine(i,seq)) || (PicardUtils.isAdenine(i,ref) && PicardUtils.isGuanine(i,seq))){	
									tempInt[0] = 1;
									//tempInt[1] = readsDepth.get(onRefCoord);
									allelePosition.put(strandRefCoord,tempInt);
								}
							}
						}
						if (refi == '-')
						{
							// It's a deletion in reference, don't advance
						}
						else
						{
							int inc = (negStrand) ? -1 : 1;
							onRefCoord += inc;
						}
						//System.err.println(onRefCoord);
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
		
		
		
		protected TreeMap<Integer,List<Byte>> alleleReads( SAMFileReader inputSam, TreeMap<Integer,Integer[]> allelePosition, boolean canPurge, String chr)
		throws Exception{
			CloseableIterator<SAMRecord> chrItForAllele = inputSam.query(chr, 0, 0, false);
			int lastBaseSeen = 0;
			int lastPurge = 0;
			int recCounter = 0;
			//List<Integer> readsDepth = new ArrayList<Integer>();
			//int j = 0;
			TreeMap<Integer,List<Byte>> alleleReadsMem = new TreeMap<Integer,List<Byte>>();
			
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
				if (unmapped || (mapQual < minMapQ))
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
					//System.err.println(seq + "\n" + ref);

					boolean negStrand = samRecord.getReadNegativeStrandFlag();
					int alignmentS = samRecord.getAlignmentStart();
					int	onRefCoord = (negStrand) ? samRecord.getUnclippedEnd() : alignmentS; 
					//int readsStart = samRecord.getUnclippedEnd() <= alignmentS ? samRecord.getUnclippedEnd() : alignmentS;
					//int readsEnd = samRecord.getUnclippedEnd() <= alignmentS ? alignmentS : samRecord.getUnclippedEnd();

					if (alignmentS < lastBaseSeen)
					{
						System.err.printf("BAM must be ordered in order: %d<%d\n",alignmentS, lastBaseSeen);
						System.exit(1);
					}
					lastBaseSeen = alignmentS;

					int seqLen = Math.min(seq.length(), ref.length());
					
					for (int i = 0; i < seqLen; i++){
						char refi = ref.charAt(i);
						//char seqi = seq.charAt(i);
						byte baseQS = (negStrand) ? baseQual[seqLen-1-i] : baseQual[i];
						Integer strandRefCoord = (negStrand) ? 0-onRefCoord : onRefCoord;
							if(allelePosition.containsKey(strandRefCoord)){
										 List<Byte> tempQuality = new ArrayList<Byte>();
										 if(alleleReadsMem.containsKey(strandRefCoord))
											 tempQuality = alleleReadsMem.get(strandRefCoord);
										tempQuality.add(baseQS);
										alleleReadsMem.put(strandRefCoord, tempQuality);
										
							}
								
							//}
						//}
						//if( onRefCoord == 2718538 && PicardUtils.isAdenine(i,seq) || PicardUtils.isGuanine(i,seq))
						//	System.err.println(baseQual[i]);
						
						if (refi == '-')
						{
							// It's a deletion in reference, don't advance
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
			return alleleReadsMem;
			
		}
		
		protected void outputSNP(TreeMap<Integer,Integer[]> allelePosition, String tableName)
		throws FileNotFoundException{
			Iterator<Integer> it = allelePosition.keySet().iterator();
			System.err.println("-----------------------------------------");
			String fn = tableName + "_SNP_all_afterBaseQfilter" + ".txt";
			PrintWriter writer = new PrintWriter(new File(fn));
			while(it.hasNext()){	
				Integer snp = it.next();
				Integer[] tempInt = allelePosition.get(snp);
				 double tempDouble1 = tempInt[0];
				 double tempDouble2 = tempInt[1];
				 double freq1 = tempDouble1/tempDouble2;
				 double freq2 = (tempDouble2-tempDouble1)/tempDouble2;
				 //System.err.println(tempDouble1);
				 //System.err.println(tempDouble2);
				 //System.err.println(freq);
				 if(tempInt[1] <= 10 && (tempInt[0] < minAlleleCount || (tempInt[1] - tempInt[0]) < minAlleleCount) ){
					 //allelePosition.remove(alleleChromPos);
					 continue;
				 }
				 else if(tempInt[1] > 10 && ((freq1 < minAlleleFreq || freq2 < minAlleleFreq) || (tempInt[0] < minAlleleCount || (tempInt[1] - tempInt[0]) < minAlleleCount))){
					 continue;
				 }
				 else{
					 writer.printf("%d\t%d\t%d\n",snp,tempInt[1],tempInt[0]);
					 
				 }
				
			}
			writer.close();
		}
		
		protected void outputSNP(TreeMap<Integer,Integer[]> allelePosition, TreeMap<Integer,List<Byte>> alleleReadsMem, String tableName)
		throws FileNotFoundException{
			Iterator<Integer> it = allelePosition.keySet().iterator();
			System.err.println("-----------------------------------------");
			String fn = tableName + "_SNP" + ".txt";
			String fnBase = tableName + "_allBaseQuality" + ".txt";
			String fnQuantBase = tableName + "_baseQuality" + ".txt";
			PrintWriter writer = new PrintWriter(new File(fn));
			PrintWriter writerQuantBase = new PrintWriter(new File(fnQuantBase));
			PrintWriter writerBase = new PrintWriter(new File(fnBase));
			while(it.hasNext()){	
				Integer snp = it.next();
				Integer[] tempInt = allelePosition.get(snp);
				List<Byte> tempByte = alleleReadsMem.get(snp);
				 double tempDouble1 = tempInt[0];
				 double tempDouble2 = tempInt[1];
				 double freq1 = tempDouble1/tempDouble2;
				 double freq2 = (tempDouble2-tempDouble1)/tempDouble2;
				 //System.err.println(tempDouble1);
				 //System.err.println(tempDouble2);
				 //System.err.println(freq);
				 if(tempInt[1] <= 10 && (tempInt[0] < minAlleleCount || (tempInt[1] - tempInt[0]) < minAlleleCount) ){
					 //allelePosition.remove(alleleChromPos);
					 continue;
				 }
				 else if(tempInt[1] > 10 && ((freq1 < minAlleleFreq || freq2 < minAlleleFreq) || (tempInt[0] < minAlleleCount || (tempInt[1] - tempInt[0]) < minAlleleCount))){
					 continue;
				 }
				 else{
					 writer.printf("%d\t%d\t%d\n",snp,tempInt[1],tempInt[0]);
					  qantileList(tempByte, writerQuantBase, snp); 
					 Iterator<Byte> listIt = tempByte.iterator();
					 while(listIt.hasNext()){
						 writerBase.printf("%d\t%d\n",snp,listIt.next());
					 }
					 
				 }
				
			}
			writer.close();
			writerQuantBase.close();
			writerBase.close();
		}
		
		protected static void qantileList(List<Byte> tempByte, PrintWriter writerBase, Integer snp){
			Byte[] a = new Byte[tempByte.size()];
			a = tempByte.toArray(a);
			Byte[] b = quicksort(a, 0, a.length-1);
			for(int i = 0 ; i < (b.length)/4 ; i++){
				writerBase.printf("%d\t%d\n",snp,b[i]);
			}
		}
		
		private static Byte[] quicksort(Byte[] array, int lo, int hi)
		{
		   if (hi > lo)
		   {
		      int partitionPivotIndex = lo;
		 
		      int newPivotIndex = partition(array, lo, hi, partitionPivotIndex);
		 
		      quicksort(array, lo, newPivotIndex-1);
		      quicksort(array, newPivotIndex+1, hi);
		   }
		   return  array;
		}

		private static int partition(Byte[] array, int lo, int hi, int pivotIndex)
		{
			Byte pivotValue = array[ pivotIndex ];
		 
		   swap(array, pivotIndex, hi); //send pivot item to the back
		 
		   int index = lo; //keep track of where the front ends
		 
		   for (int i = lo; i < hi; i++) //check from the front to the back
		   {
		      //swap if the current value is less than the pivot
		      if ( (array[i]).compareTo(pivotValue) <= 0 )
		      {
		         swap(array, i, index);
		         index++;
		      }
		   }
		 
		   swap(array, hi, index); //put pivot item in the middle
		 
		   return index;
		}
		
		private static void swap(Byte[] array, int i, int j)
		   {
			  Byte temp = array[i];
		      array[i] = array[j];
		      array[j] = temp;
		   }


}
