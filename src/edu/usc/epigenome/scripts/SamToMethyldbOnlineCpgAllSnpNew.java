package edu.usc.epigenome.scripts;

import net.sf.samtools.*;
import net.sf.samtools.util.CloseableIterator;


import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.sql.*;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.TreeSet;
import java.util.SortedMap;
import java.util.TreeMap;


import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;

import edu.usc.epigenome.genomeLibs.PicardUtils;
import edu.usc.epigenome.genomeLibs.MethylDb.Cytosine;



public class SamToMethyldbOnlineCpgAllSnpNew {

	/**
	 * @param args
	 */
		
		final private static String prefix = "methylCGsRich_ASM_AllSnp_";
		final private static String USAGE = "SamToMethyldbOnlineCpgAllSnpNew [opts] chr sampleName snpfile.txt file1.bam";

		final private static int PURGE_INTERVAL = 20000; // We purge our stored Cpgs once we get this many bases past them.
		//final private static int ALLELE_GA_NUMBER = 1; //for each reads, when there are more than 1 GA position is different from reference sequence, we define it belongs to another allele. 

		public static String connStr = "jdbc:mysql://epifire2.epigenome.usc.edu/asm_cr";
		//public static String connStr = "jdbc:mysql://hpc2721/asm_cr";
		//public static String connStr = null;
		protected static Connection cConn = null;
		//mysql_db_server: epifire2.epigenome.usc.edu
		
		
		/**
		 * object vars
		 */

		
		/**
		 * @param args
		 */

		@Option(name="-minConv",usage="minimum number of converted cytosines required")
		protected int minConv = 1;
		@Option(name="-useCpgsToFilter",usage=" Use CpGs and CpHs to filter if true, otherwise just CpHs (default false)")
		protected boolean useCpgsToFilter = false;
		@Option(name="-outputHcphs",usage=" Output HCpH cytosines (can't use more than 1 input file)")
		protected boolean outputHcphs = false;
		@Option(name="-minMapQ",usage="minimum mapping quality (default 30)")
		protected int minMapQ = 20;
		@Option(name="-minBaseQual",usage="minimum Base quality (default 10)")
		protected int minBaseQual = 10;
		@Option(name="-minCpgBaseQual",usage="minimum CpG Base quality (default 14)")
		protected int minCpgBaseQual = 10;
		@Option(name="-minAlleleCount",usage="minimum Allele Count (default 3)")
		protected static int minAlleleCount = 1;
		@Option(name="-minAlleleFreq",usage="minimum BAllele Frequency (default 0.2)")
		protected static double minAlleleFreq = 0.10;

		@Option(name="-debug",usage=" Debugging statements (default false)")
		protected boolean debug = false;

		private boolean asmFlag = true;
		
		// receives other command line parameters than options
		@Argument
		private List<String> stringArgs = new ArrayList<String>();

		
		
		/**
		 * @param args
		 */
		public static void main(String[] args)
		throws Exception
		{
			new SamToMethyldbOnlineCpgAllSnpNew().doMain(args);
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
				
				if (this.outputHcphs && (stringArgs.size()>2)) throw new CmdLineException("Can't use -outputHcphs with multiple input bam files");
			}
			catch (CmdLineException e)
			{
				System.err.println(e.getMessage());
				// print the list of available options
				parser.printUsage(System.err);
				System.err.println();
				return;
			}
			String chr = stringArgs.remove(0);
			String sampleName = stringArgs.remove(0);
			String snpFileName = stringArgs.remove(0);
			String fn = stringArgs.remove(0);
			BufferedReader br = new BufferedReader(new FileReader(snpFileName));
			TreeSet<Integer> snpPosition = new TreeSet<Integer>();
			String line;
			//int a = 0;
			while( (line = br.readLine()) != null){
				String[] tmpArray = line.split("\t");
				int tmpSnp = Integer.parseInt(tmpArray[0]);
				if(snpPosition.contains(tmpSnp-1)){
					snpPosition.remove(tmpSnp-1);
					continue;
				}
				if(snpPosition.contains(tmpSnp-2)){
					snpPosition.remove(tmpSnp-2);
					continue;
				}
				if(snpPosition.contains(tmpSnp-3)){
					snpPosition.remove(tmpSnp-3);
					continue;
				}
				if(snpPosition.contains(tmpSnp-4)){
					snpPosition.remove(tmpSnp-4);
					continue;
				}
				if(snpPosition.contains(tmpSnp-5)){
					snpPosition.remove(tmpSnp-5);
					continue;
				}
				snpPosition.add(tmpSnp);
				//System.err.println(Integer.parseInt(tmpArray[0]));
				//a++;
			}
			br.close();
			Iterator<Integer> it = snpPosition.iterator();
			if(connStr == null){
				BufferedReader hostBr = new BufferedReader(new FileReader("/home/uec-00/shared/production/database/mysql.host.txt"));
				String hostName = hostBr.readLine();
				connStr = "jdbc:mysql://" + hostName + "/asm_cr";
			}
			
			setupDb(connStr);
			int recCounter = 0;
			SortedMap<String,Cytosine> cytosines = new TreeMap<String,Cytosine>();
			String tableName = prefix + sampleName + "_" + chr;
			Cytosine.outputChromToDb(cytosines, tableName, cConn, asmFlag);
			File inputSamOrBamFile = new File(fn);
			
			
			while(it.hasNext()){
				int snp = it.next();
				SAMFileReader inputSam = new SAMFileReader(inputSamOrBamFile);
				inputSam.setValidationStringency(SAMFileReader.ValidationStringency.SILENT);
				CloseableIterator<SAMRecord> chrIt = inputSam.queryOverlapping(chr, snp, snp);
				char A_BaseUpperCase = '0';
				char B_BaseUpperCase = '0';
				boolean agFlag = false;
				boolean agNegFlag = false;
				while (chrIt.hasNext()){
					SAMRecord samRecord = chrIt.next();
					int mapQual = samRecord.getMappingQuality();
					byte[] baseQual = samRecord.getBaseQualities();
					boolean unmapped = samRecord.getReadUnmappedFlag();
					if (unmapped || (mapQual <= minMapQ))
					{
						continue;
					}
					String seq = PicardUtils.getReadString(samRecord, true);
					try
					{
						String ref = PicardUtils.refStr(samRecord, true);
						boolean negStrand = samRecord.getReadNegativeStrandFlag();
						int alignmentS = samRecord.getAlignmentStart();
						int	onRefCoord = (negStrand) ? samRecord.getUnclippedEnd() : alignmentS; 
						int pos = (negStrand) ? onRefCoord - snp : snp - onRefCoord;
						int seqLen = Math.min(seq.length(), ref.length());
						byte baseQS = (negStrand) ? baseQual[seqLen-1-pos] : baseQual[pos];
						if(baseQS > minBaseQual){
							if((PicardUtils.isGuanine(pos,ref) && PicardUtils.isAdenine(pos,seq)) || (PicardUtils.isAdenine(pos,ref) && PicardUtils.isGuanine(pos,seq))){
								agFlag = true;

									A_BaseUpperCase = ref.charAt(pos);
								
							
									B_BaseUpperCase = seq.charAt(pos);
								
								if(negStrand){
									agNegFlag = true;
								}
								
							}
							else if ((PicardUtils.isGuanine(pos,ref) && PicardUtils.isCytosine(pos, seq, true)) || (PicardUtils.isAdenine(pos,ref) && PicardUtils.isCytosine(pos, seq, true)) || (PicardUtils.isGuanine(pos,seq) && PicardUtils.isCytosine(pos, ref, true)) || (PicardUtils.isAdenine(pos,seq) && PicardUtils.isCytosine(pos, ref, true))){	
								//alleleNum++;// allele reads number
								if(A_BaseUpperCase == '0'){
									A_BaseUpperCase = ref.charAt(pos);
								}
								if(B_BaseUpperCase == '0'){
									B_BaseUpperCase = seq.charAt(pos);
								}
							}
						}
					}
					catch (Exception e)
					{
						System.err.println("-----------------------------------------");
						System.err.println("Couldn't handle seq #" + recCounter);
						System.err.println(seq);
						e.printStackTrace(System.err);
						System.err.println("-----------------------------------------");
					}
					
				}	
				chrIt.close();
				chrIt = inputSam.queryOverlapping(chr, snp, snp);
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
					recCounter++;				
					if ((recCounter % 1E5)==0)
					{
						System.err.printf("On new record #%d\n",recCounter);
						 
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
						int numConverted = 0;
						int convStart = Integer.MAX_VALUE;
						int seqLen = Math.min(seq.length(), ref.length());
						
						int pos = (negStrand) ? onRefCoord - snp : snp - onRefCoord;
						if(pos < 0 ){
							System.err.println("snp position error!");
						}
						boolean seqFlag = true;
						 if( (PicardUtils.isGuanine(pos,ref) && PicardUtils.isAdenine(pos,seq)) || (PicardUtils.isAdenine(pos,ref) && PicardUtils.isGuanine(pos,seq)) || (PicardUtils.isGuanine(pos,ref) && PicardUtils.isCytosine(pos, seq, true)) || (PicardUtils.isAdenine(pos,ref) && PicardUtils.isCytosine(pos, seq, true)) || (PicardUtils.isGuanine(pos,seq) && PicardUtils.isCytosine(pos, ref, true)) || (PicardUtils.isAdenine(pos,seq) && PicardUtils.isCytosine(pos, ref, true))){
							 seqFlag = false;
						 }
						 if(agFlag){
							 if(agNegFlag && !negStrand){
								 continue record;
							 }
							 else if(!agNegFlag && negStrand){
								 continue record;
							 }
						 }
						//boolean iscpgOverlapSnp = PicardUtils.isCpg(pos,ref);
						//boolean iscpgOverlapSnpAdj = false;
						//if(pos > 0 && pos < seqLen-1){
						//	iscpgOverlapSnpAdj = (negStrand) ? PicardUtils.isCpg(pos+1,ref) : PicardUtils.isCpg(pos-1,ref);
						//}
						
						for (int i = 0; i < seqLen; i++)
						{
							
							char refi = ref.charAt(i);
							char seqi = seq.charAt(i);
							char preBaseRef = PicardUtils.preBaseRef(i, ref);
							char preBaseSeq = PicardUtils.preBaseSeq(i, seq);
							char nextBaseRef = PicardUtils.nextBaseRef(i, ref);
							char nextBaseSeq = PicardUtils.nextBaseSeq(i, seq);
							
							if (PicardUtils.isCytosine(i,ref,false) && PicardUtils.isCytosine(i,seq,true)) // The last one is too tricky to deal with since we don't know context
							{
								
								boolean ishch = PicardUtils.isHch(i,ref);
								boolean iscpg = PicardUtils.isCpg(i,ref);
								
								boolean conv = PicardUtils.isConverted(i,ref,seq);
								
								byte baseCpgQS = (negStrand) ? baseQual[seqLen-1-i] : baseQual[i];
								byte baseQS = (negStrand) ? baseQual[seqLen-1-pos] : baseQual[pos];
								
								if (conv && (this.useCpgsToFilter || (!iscpg)) && baseCpgQS > minCpgBaseQual) numConverted++;

								if ((convStart==Integer.MAX_VALUE) && (numConverted>=this.minConv) && baseCpgQS > minCpgBaseQual)
								{
									convStart = i;
								}
								
								if(!outputHcphs & ishch){
									
								}
								else if( iscpg && baseCpgQS > minCpgBaseQual && baseQS > minBaseQual ){
									if(pos - i == 1 || pos - i == 0){
										continue;
									}
									Cytosine cytosine = findOrCreateCytosine(cytosines, onRefCoord, negStrand, preBaseRef, nextBaseRef, A_BaseUpperCase, B_BaseUpperCase, snp);
									if (cytosine.getNextBaseRef() == '0') cytosine.setNextBaseRef(nextBaseRef);
									if (cytosine.getPreBaseRef() == '0') cytosine.setPreBaseRef(preBaseRef);
									if (cytosine.getA_BaseUpperCase() == '0') cytosine.setA_BaseUpperCase(A_BaseUpperCase);
									if (cytosine.getB_BaseUpperCase() == '0') cytosine.setB_BaseUpperCase(B_BaseUpperCase);
									this.incrementCytosine(cytosine, seqi, i<convStart, preBaseSeq, nextBaseSeq, seqFlag);
								}
							}
							
							
							boolean oppositeCpg = PicardUtils.isOppositeCpg(i,ref);
							if (oppositeCpg)
							{
								// Look for cpg on opposite strand
								// Be careful, the "nextBaseRef" is now on the opposite strand!!
								char nextBaseRefRev = PicardUtils.nextBaseRef(i, ref, true);
								char preBaseRefRev = PicardUtils.preBaseRef(i, ref, true);
								byte baseCpaQS = (negStrand) ? baseQual[seqLen-1-i] : baseQual[i];
								if(baseCpaQS > minBaseQual){
									 Cytosine cytosine = findOrCreateCytosine(cytosines, onRefCoord, !negStrand, preBaseRefRev, nextBaseRefRev, A_BaseUpperCase, B_BaseUpperCase, snp);
									 this.incrementOppositeCpg(cytosine, seqi);
									
								}
							}
							

							// Increment genomic coord position
							if (refi == '-')
							{
								// It's a deletion in reference, don't advance
							}
							else
							{
								int inc = (negStrand) ? -1 : 1;
								onRefCoord += inc;
							}

						} // i (pos within read)
						
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

				chrIt.close();
				inputSam.close();
				
			}
				// And output the chrom
			if (cytosines != null) 	Cytosine.outputCytocinesToDb(cytosines, tableName, asmFlag);

			cleanupDb();
			
		}


		protected static Cytosine findOrCreateCytosine(Map<String,Cytosine> cytosines, int onRefCoord, boolean negStrand, char preBaseRef, char nextBaseRef, char A_BaseUpperCase, char B_BaseUpperCase, int alleleChromPos)
		{
			
			String index = Integer.toString(onRefCoord) + "~" + Integer.toString(alleleChromPos);
			Cytosine cytosine = cytosines.get(index);
			
			if (cytosine == null)
			{
				cytosine = new Cytosine(onRefCoord,negStrand,alleleChromPos);
				cytosine.setNextBaseRef(nextBaseRef);
				cytosine.setPreBaseRef(preBaseRef);
				cytosine.setA_BaseUpperCase(A_BaseUpperCase);
				cytosine.setB_BaseUpperCase(B_BaseUpperCase);
				
				cytosines.put(index, cytosine);
			}
			return cytosine;
			
		}
		
		protected void incrementCytosine(Cytosine cytosine, char seqChar, boolean nonconvFilter, char preBaseSeq, char nextBaseSeq, boolean seqFlag) 
		throws Exception
		{
			int totalReads = 0, cReads = 0, tReads = 0, cReadsNonconvFilt = 0, agReads = 0, preBaseGreads = 0, preBaseTotalReads = 0, nextBaseGreads = 0, nextBaseTotalReads = 0, A_CReads = 0, B_CReads = 0, A_TReads = 0, B_TReads = 0; 
			
			switch (seqChar)
			{
			case 'N':
			case '0':
				break;
			case 'A':
			case 'G':
				agReads = 1;
				totalReads = 1;
				break;
			case 'T':
				tReads = 1;
				totalReads = 1;
				if ( seqFlag ) A_TReads = 1; 
				else B_TReads = 1;
				break;
			case 'C':
				if (nonconvFilter) cReadsNonconvFilt = 1; else cReads = 1;
				totalReads = 1;
				if ( seqFlag ) A_CReads = 1;
				else B_CReads = 1;
				break;
			default:
				throw new Exception("Can't recognize seq char: " + seqChar);
			}
		
			// And the previous base
			switch (preBaseSeq)
			{
			case 'N':
			case '0':
				break;
			case 'G':
				preBaseTotalReads = 1;
				preBaseGreads = 1;
				break;
			case 'T':
			case 'A':
			case 'C':
				preBaseTotalReads = 1;
				break;
			default:
				throw new Exception("Can't recognize seq char: " + seqChar);
			}
			
			// And the next base
			switch (nextBaseSeq)
			{
			case 'N':
			case '0':
				break;
			case 'G':
				nextBaseTotalReads = 1;
				nextBaseGreads = 1;
				break;
			case 'T':
			case 'A':
			case 'C':
				nextBaseTotalReads = 1;
				break;
			default:
				throw new Exception("Can't recognize seq char: " + seqChar);
			}
		
			cytosine.totalReads += totalReads;
			cytosine.cReads += cReads;
			cytosine.tReads += tReads;
			cytosine.cReadsNonconversionFilt += cReadsNonconvFilt;
			cytosine.agReads += agReads;
			
			cytosine.preBaseGreads += preBaseGreads;
			cytosine.preBaseTotalReads += preBaseTotalReads;
			cytosine.nextBaseGreads += nextBaseGreads;
			cytosine.nextBaseTotalReads += nextBaseTotalReads;
			
			cytosine.A_CReads += A_CReads; 
			cytosine.B_CReads +=B_CReads;
			cytosine.A_TReads += A_TReads;
			cytosine.B_TReads += B_TReads;
		}
		
		protected void incrementOppositeCpg(Cytosine cytosine, char seqChar) 
		throws Exception
		{
			int aReadsOpposite = 0, totalReadsOpposite = 0;
			
			
			switch (seqChar)
			{
			case 'N':
			case '0':
				break;
			case 'A':
				aReadsOpposite = 1;
				totalReadsOpposite = 1;
				break;
			case 'G':
			case 'T':
			case 'C':
				totalReadsOpposite = 1;
				break;
			default:
				throw new Exception("Can't recognize seq char: " + seqChar);
			}

			cytosine.aReadsOpposite += aReadsOpposite;
			cytosine.totalReadsOpposite += totalReadsOpposite;
		}
		

		protected static void setupDb(String connStrCytosine)
		throws Exception
		{
			if (cConn == null)
			{
				String connStr = connStrCytosine;
				//String connStr = MethylDbQuerier.connStr;
				Class.forName("com.mysql.jdbc.Driver").newInstance();
				System.err.println("Getting connection for " + connStr);
				cConn = DriverManager.getConnection(connStrCytosine, "yaping", "lyping1986");
			}
			
		}
		
		protected static void cleanupDb()
		throws Exception
		{
			cConn.close();
		}
}
