package edu.usc.epigenome.scripts;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.PrintWriter;
import java.sql.*;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.SortedMap;
import java.util.TreeMap;
import java.util.logging.Logger;

import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.util.CloseableIterator;

import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;

import edu.usc.epigenome.genomeLibs.PicardUtils;
import edu.usc.epigenome.genomeLibs.MethylDb.Cytosine;
import edu.usc.epigenome.genomeLibs.MethylDb.MethylDbUtils;

public class SamToMethyldbOnlineAllCytocine {

	/**
	 * @param args
	 */
		
		final private static String prefix = "methylCGsRich_";
		final private static String USAGE = "SamToMethyldbOnline [opts] sampleName file1.bam file2.bam ...";

		final private static int PURGE_INTERVAL = 20000; // We purge our stored Cpgs once we get this many bases past them.
		public static String connStr = "jdbc:mysql://epifire2.epigenome.usc.edu/gnome_seq";
		//public static String connStr = "jdbc:mysql://hpc2721/gnome_seq";
		//public static String connStr = null;
		protected static Connection cConn = null;
		//mysql_db_server: epifire2.epigenome.usc.edu
		
		/**
		 * object vars
		 */

		
		/**
		 * @param args
		 */
		@Option(name="-chrom",multiValued=true,usage="One or more chroms, eg. --chrom chr1 --chrom chr5")
		protected List<String> chrs = new ArrayList<String>(25);
		@Option(name="-minConv",usage="minimum number of converted cytosines required")
		protected int minConv = 1;
//		@Option(name="-numCycles",usage="Number of cycles to track")
//		protected int numCycles = 100;
//		@Option(name="-outputReads",usage=" Outputs one line per read (default false)")
//		protected boolean outputReads = false;
		@Option(name="-useCpgsToFilter",usage=" Use CpGs and CpHs to filter if true, otherwise just CpHs (default false)")
		protected boolean useCpgsToFilter = false;
		@Option(name="-useGchsToFilter",usage=" Use GCHs and HCpHs to filter if true, otherwise just HCpHs (default false)")
		protected boolean useGchsToFilter = false;
		@Option(name="-useHcgsToFilter",usage=" Use HCGs and HCpHs to filter if true, otherwise just HCpHs (default false)")
		protected boolean useHcgsToFilter = false;
		@Option(name="-outputHcphs",usage=" Output HCpH cytosines (can't use more than 1 input file)")
		protected boolean outputHcphs = false;
		@Option(name="-minMapQ",usage="minimum mapping quality (default 30)")
		protected int minMapQ = 30;
		@Option(name="-minHcphCoverage",usage="the minimum number of total reads to include a HCph (default 10)")
		protected int minHcphCoverage = 10;
		@Option(name="-minHcphFrac",usage="minimum methylation fraction to include a HCph")
		protected double minHcphFrac = 0.2;
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
			new SamToMethyldbOnlineAllCytocine().doMain(args);
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

			String sampleName = stringArgs.remove(0);

			int recCounter = 0;
			int usedCounter = 0;
			int filteredOutCounter = 0;
			int gchFilteredOutCounter = 0;
			int gchUsedCounter = 0;
			int hcgFilteredOutCounter = 0;
			int hcgUsedCounter = 0;
			// Cph counters
			int numCphConvertedWithFilt = 0;
			int numCphTotalWithFilt = 0;
			int numCphConvertedNoFilt = 0;
			int numCphTotalNoFilt = 0;
			// Hcph counters
			int numHcphConvertedWithFilt = 0;
			int numHcphTotalWithFilt = 0;
			int numHcphConvertedNoFilt = 0;
			int numHcphTotalNoFilt = 0;
			//CpG counters
			int numCpgConvertedWithFilt = 0;
			int numCpgTotalWithFilt = 0;
			int numCpgConvertedNoFilt = 0;
			int numCpgTotalNoFilt = 0;
			//Gch counters
			int numGchConvertedWithFilt = 0;
			int numGchTotalWithFilt = 0;
			int numGchConvertedNoFilt = 0;
			int numGchTotalNoFilt = 0;
			//Hcg counters
			int numHcgConvertedWithFilt = 0;
			int numHcgTotalWithFilt = 0;
			int numHcgConvertedNoFilt = 0;
			int numHcgTotalNoFilt = 0;
			
			if(connStr == null){
				BufferedReader hostBr = new BufferedReader(new FileReader("/home/uec-00/shared/production/database/mysql.host.txt"));
				String hostName = hostBr.readLine();
				connStr = "jdbc:mysql://" + hostName + "/gnome_seq";
			}
			
			setupDb(connStr);
			
			// Iterate through chroms
			if (chrs.size()==0) chrs = MethylDbUtils.CHROMS;
			for (final String chr : chrs)
			{
				SortedMap<Integer,Cytosine> cytosines = new TreeMap<Integer,Cytosine>();
				String tableName = prefix + sampleName + "_" + chr;
				Cytosine.outputChromToDb(cytosines, tableName, cConn, this.minHcphCoverage, this.minHcphFrac);

				
				for (final String fn : stringArgs)
				{
					File inputSamOrBamFile = new File(fn);

					final SAMFileReader inputSam = new SAMFileReader(inputSamOrBamFile);
					inputSam.setValidationStringency(SAMFileReader.ValidationStringency.SILENT);
					
					CloseableIterator<SAMRecord> chrIt = inputSam.query(chr, 0, 0, false);
					
					// We can only purge if we are the only input file and we are sorted.
					boolean canPurge = ((inputSam.hasIndex()) && (stringArgs.size() == 1));
					Logger.getLogger(Logger.GLOBAL_LOGGER_NAME).info(String.format("Able to purge at interval %d? %s\n", 
							PURGE_INTERVAL, (canPurge) ? "Yes" : "Purging not available - either unsorted BAM or multiple BAMs for the same chrom"));
					int lastBaseSeen = 0;
					int lastPurge = 0;
					record: while (chrIt.hasNext())
					{
						if (canPurge && (lastBaseSeen > (lastPurge+PURGE_INTERVAL)))
						{
							//System.err.printf("On base %d, purging everything before %d\n", lastBaseSeen,lastPurge);
							Cytosine.outputCytocinesToDb(cytosines.headMap(new Integer(lastPurge)), tableName, this.minHcphCoverage, this.minHcphFrac);
							
							// Weird, if i just set cytosines to be the tailMap (as in 1, below) garbage collection doesn't actually clean up
							// the old part (just backed by the original data structure). So i actually have to copy it to a new map. (as in 2)
							//(1) cytosines = cytosines.tailMap(new Integer(lastPurge));
							cytosines = new TreeMap<Integer,Cytosine>(cytosines.tailMap(new Integer(lastPurge))); // (2)
							
							lastPurge = lastBaseSeen;
						}
						
						SAMRecord samRecord = chrIt.next();

						// Filter low qual
						int mapQual = samRecord.getMappingQuality();
						boolean unmapped = samRecord.getReadUnmappedFlag();
						if (unmapped || (mapQual < minMapQ))
						{
							continue record;
						}


						String seq = PicardUtils.getReadString(samRecord, true);

						recCounter++;
						if ((recCounter % 1E5)==0)
						{
							System.err.printf("On new record #%d, purged tree size:%d\n",recCounter,cytosines.size());
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

							if ((this.outputHcphs) && (alignmentS < lastBaseSeen))
							{
								System.err.printf("BAM must be ordered in order to handle -outputCphs: %d<%d\n",alignmentS, lastBaseSeen);
								System.exit(1);
							}
							lastBaseSeen = alignmentS;
								
								

							int numConverted = 0;
							int convStart = Integer.MAX_VALUE;
							int gchNumConverted = 0;
							int gchConvStart = Integer.MAX_VALUE;
							int hcgNumConverted = 0;
							int hcgConvStart = Integer.MAX_VALUE;
							
							int seqLen = Math.min(seq.length(), ref.length());
							for (int i = 0; i < seqLen; i++)
							{
								char refi = ref.charAt(i);
								char seqi = seq.charAt(i);
								char preBaseRef = PicardUtils.preBaseRef(i, ref);
								char preBaseSeq = PicardUtils.preBaseSeq(i, seq);
								char nextBaseRef = PicardUtils.nextBaseRef(i, ref);
								char nextBaseSeq = PicardUtils.nextBaseSeq(i, seq);

								//if ((i < (seqLen-1)) && PicardUtils.isCytosine(i,ref)) // The last one is too tricky to deal with since we don't know context
								if (PicardUtils.isCytosine(i,ref,false) && PicardUtils.isCytosine(i,seq,true)) // The last one is too tricky to deal with since we don't know context
								{
									boolean isgch = PicardUtils.isGch(i,ref);
									boolean isgcg = PicardUtils.isGcg(i,ref);
									boolean ishcg = PicardUtils.isHcg(i,ref);
									boolean ishch = PicardUtils.isHch(i,ref);
									boolean iscpg = PicardUtils.isCpg(i,ref);
									
									boolean conv = PicardUtils.isConverted(i,ref,seq);
									
//									if (iscpg && (nextBaseRef != 'G'))
//									{
//										System.err.printf("Cpg status and next base don't match:\nref=%s\nseq=%s\nref_%d=%c, seq_%d=%c\n",
//												ref, seq, i, refi, i, seqi);
//									}
									

									//if (conv && this.useCpgsToFilter || !iscpg) numConverted++;
									if (conv && (this.useCpgsToFilter || (!iscpg && (i < (seqLen-1)) && nextBaseSeq != 'G'))) numConverted++;
									if (conv && (this.useGchsToFilter || (ishch && preBaseSeq != 'G' && nextBaseSeq != 'G'))) gchNumConverted++;
									if (conv && (this.useHcgsToFilter || (!ishcg && (preBaseSeq == 'G' || nextBaseSeq != 'G')))) hcgNumConverted++;
									// If this is the first legal one , note it
									if ((convStart==Integer.MAX_VALUE) && (numConverted>=this.minConv) && (i < (seqLen-1)) )
									{
										convStart = i;
									}
									if ((gchConvStart==Integer.MAX_VALUE) && (gchNumConverted>=this.minConv) )
									{
										gchConvStart = i;
									}
									if ((hcgConvStart==Integer.MAX_VALUE) && (hcgNumConverted>=this.minConv) )
									{
										hcgConvStart = i;
									}
									
									if(!outputHcphs && ishch){
										
									}
									else{
										Cytosine cytosine = findOrCreateCytosine(cytosines, onRefCoord, negStrand, preBaseRef, nextBaseRef);
										if (cytosine.getNextBaseRef() == '0') cytosine.setNextBaseRef(nextBaseRef);
										if (cytosine.getPreBaseRef() == '0') cytosine.setPreBaseRef(preBaseRef);
										this.incrementCytosine(cytosine, seqi, i<convStart || i<gchConvStart, preBaseSeq, nextBaseSeq);
									}
									
									

									//if(iscpg)
									if (iscpg && (i < (seqLen-1)) && nextBaseSeq == 'G')
									{
										if (i<convStart)
										{
											// In the non-conversion filter zone
											filteredOutCounter++;
											//System.err.printf("Rec %d\tpos=%d\n",recCounter,i);
											numCpgTotalNoFilt++;
											if (conv) numCpgConvertedNoFilt++;
										}
										else
										{
											// Past the non-conversion filter, use it
											usedCounter++;
											numCpgTotalWithFilt++;
											if (conv) numCpgConvertedWithFilt++;
										}

									}
									
									//if(isgch)
									if (isgch && preBaseSeq == 'G' && nextBaseSeq != 'G')
									{
										if (i<gchConvStart)
										{
											// In the non-conversion filter zone
											gchFilteredOutCounter++;
											//System.err.printf("Rec %d\tpos=%d\n",recCounter,i);
											numGchTotalNoFilt++;
											if (conv) numGchConvertedNoFilt++;
										}
										else
										{
											// Past the non-conversion filter, use it
											gchUsedCounter++;
											numGchTotalWithFilt++;
											if (conv) numGchConvertedWithFilt++;
										}

									}
									
									//if(ishcg)
									if (ishcg && preBaseSeq != 'G' && nextBaseSeq == 'G')
									{
										if (i<hcgConvStart)
										{
											// In the non-conversion filter zone
											hcgFilteredOutCounter++;
											//System.err.printf("Rec %d\tpos=%d\n",recCounter,i);
											numHcgTotalNoFilt++;
											if (conv) numHcgConvertedNoFilt++;
										}
										else
										{
											// Past the non-conversion filter, use it
											hcgUsedCounter++;
											numHcgTotalWithFilt++;
											if (conv) numHcgConvertedWithFilt++;
										}

									}
									
									if (!iscpg && (i < (seqLen-1)) && nextBaseSeq != 'G'){
										
										numCphTotalNoFilt++;
										if (conv) numCphConvertedNoFilt++;
										
										if (i>=convStart)
										{
											numCphTotalWithFilt++;
											if (conv) numCphConvertedWithFilt++;
										}
									}
									
									if (ishch && preBaseSeq != 'G' && nextBaseSeq != 'G'){
										
										numHcphTotalNoFilt++;
										if (conv) numHcphConvertedNoFilt++;
										
										if (i>=gchConvStart)
										{
											numHcphTotalWithFilt++;
											if (conv) numHcphConvertedWithFilt++;
										}
									}



								} // IsCytosine

								boolean oppositeGch = PicardUtils.isOppositeGch(i,ref);
								boolean oppositeGcg = PicardUtils.isOppositeGcg(i,ref);
								boolean oppositeHcg = PicardUtils.isOppositeHcg(i,ref);
								boolean oppositeHch = PicardUtils.isOppositeHch(i,ref);
								
								//if( oppositeGcg || oppositeHcg )
								if (oppositeGch || oppositeGcg || oppositeHcg || (oppositeHch & outputHcphs))
								{
		//							// Look for cpg on opposite strand
			//						// Be careful, the "nextBaseRef" is now on the opposite strand!!
									
									char nextBaseRefRev = PicardUtils.nextBaseRef(i, ref, true);
									char preBaseRefRev = PicardUtils.preBaseRef(i, ref, true);
									Cytosine oppositeCytosine = findOrCreateCytosine(cytosines, onRefCoord, !negStrand, preBaseRefRev, nextBaseRefRev);
							//		
									this.incrementOppositeCytosine(oppositeCytosine, seqi);
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
//							chrIt.close();
//							System.exit(1);
						}

					} // record

					chrIt.close();
					inputSam.close();

				}

				// And output the chrom
				if (cytosines != null) 	Cytosine.outputCytocinesToDb(cytosines, tableName, this.minHcphCoverage, this.minHcphFrac);

			}
			
			
			double frac = (double)filteredOutCounter/((double)usedCounter+(double)filteredOutCounter);
			double gchFrac = (double)gchFilteredOutCounter/((double)gchUsedCounter+(double)gchFilteredOutCounter);
			double hcgFrac = (double)hcgFilteredOutCounter/((double)hcgUsedCounter+(double)hcgFilteredOutCounter);
			System.err.printf("Lost %f%% due to CpG non-converion filter\n%d CpGs filtered for non-conversion, %d CpGs used (MinConv=%d,UseCpgs=%s)\n",
					frac*100.0, filteredOutCounter, usedCounter, this.minConv, String.valueOf(this.useCpgsToFilter));
			System.err.printf("Lost %f%% due to Hcg non-converion filter\n%d Hcgs filtered for non-conversion, %d Hcgs used (MinConv=%d,UseHcgs=%s)\n",
					hcgFrac*100.0, hcgFilteredOutCounter, hcgUsedCounter, this.minConv, String.valueOf(this.useHcgsToFilter));
			System.err.printf("Lost %f%% due to Gch non-converion filter\n%d Gchs filtered for non-conversion, %d Gchs used (MinConv=%d,UseGchs=%s)\n",
					gchFrac*100.0, gchFilteredOutCounter, gchUsedCounter, this.minConv, String.valueOf(this.useGchsToFilter));
			System.err.printf("Found %d reads total\n", recCounter);
			System.err.printf("CpH conversion rate: before filter=%f, after filter=%f\n",
					100.0 * ((double)numCphConvertedNoFilt/(double)numCphTotalNoFilt),
					100.0 * ((double)numCphConvertedWithFilt/(double)numCphTotalWithFilt));
			System.err.printf("numCphConvertedNoFilt: %d, numCphTotalNoFilt: %d\n",
					numCphConvertedNoFilt, numCphTotalNoFilt);
			System.err.printf("HCpH conversion rate: before filter=%f, after filter=%f\n",
					100.0 * ((double)numHcphConvertedNoFilt/(double)numHcphTotalNoFilt),
					100.0 * ((double)numHcphConvertedWithFilt/(double)numHcphTotalWithFilt));
			System.err.printf("CpG conversion rate: before filter=%f, after filter=%f\n",
					100.0 * ((double)numCpgConvertedNoFilt/(double)numCpgTotalNoFilt),
					100.0 * ((double)numCpgConvertedWithFilt/(double)numCpgTotalWithFilt));
			System.err.printf("HCG conversion rate: before filter=%f, after filter=%f\n",
					100.0 * ((double)numHcgConvertedNoFilt/(double)numHcgTotalNoFilt),
					100.0 * ((double)numHcgConvertedWithFilt/(double)numHcgTotalWithFilt));
			System.err.printf("GCH conversion rate: before filter=%f, after filter=%f\n",
					100.0 * ((double)numGchConvertedNoFilt/(double)numGchTotalNoFilt),
					100.0 * ((double)numGchConvertedWithFilt/(double)numGchTotalWithFilt));
		}


		protected static Cytosine findOrCreateCytosine(Map<Integer,Cytosine> cytosines, int onRefCoord, boolean negStrand, char preBaseRef, char nextBaseRef)
		{
			Cytosine cytosine = cytosines.get(new Integer(onRefCoord));
			if (cytosine == null)
			{
				cytosine = new Cytosine(onRefCoord,negStrand);
				cytosine.setNextBaseRef(nextBaseRef);
				cytosine.setPreBaseRef(preBaseRef);
				
				cytosines.put(new Integer(onRefCoord), cytosine);
			}
			return cytosine;
		}
		
		
		
		protected void incrementOppositeCytosine(Cytosine cytosine, char seqChar) 
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
		
		protected void incrementCytosine(Cytosine cytosine, char seqChar, boolean nonconvFilter, char preBaseSeq, char nextBaseSeq) 
		throws Exception
		{
			int totalReads = 0, cReads = 0, tReads = 0, cReadsNonconvFilt = 0, agReads = 0, preBaseGreads = 0, preBaseTotalReads = 0, nextBaseGreads = 0, nextBaseTotalReads = 0;
			
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
				break;
			case 'C':
				if (nonconvFilter) cReadsNonconvFilt = 1; else cReads = 1;
				totalReads = 1;
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
