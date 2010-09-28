package edu.usc.epigenome.scripts;

import net.sf.samtools.*;
import net.sf.samtools.util.CloseableIterator;


import java.io.File;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.SortedSet;
import java.util.TreeSet;
import java.util.SortedMap;
import java.util.TreeMap;
import java.util.logging.Logger;

import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;

import edu.usc.epigenome.genomeLibs.PicardUtils;
import edu.usc.epigenome.genomeLibs.MethylDb.Cytosine;
import edu.usc.epigenome.genomeLibs.MethylDb.MethylDbUtils;


public class SamToMethyldbOfflineAllCytocineASM {

	/**
	 * @param args
	 */
		
		final private static String prefix = "methylCGsRich_ASM_";
		final private static String USAGE = "SamToMethyldbOfflineAllCytocineASM [opts] sampleName file1.bam file2.bam ...";

		final private static int PURGE_INTERVAL = 20000; // We purge our stored Cpgs once we get this many bases past them.
		final private static int ALLELE_GA_NUMBER = 1; //for each reads, when there are more than 1 GA position is different from reference sequence, we define it belongs to another allele. 

		
		/**
		 * object vars
		 */

		
		/**
		 * @param args
		 */
		@Option(name="-chrom",multiValued=true,usage="One or more chroms, eg. --chrom chr1 --chrom chr5")
		protected List<String> chrs = new ArrayList<String>(25);
		@Option(name="-minConv",usage="minimum number of converted cytosines required")
		protected int minConv = 0;
//		@Option(name="-numCycles",usage="Number of cycles to track")
//		protected int numCycles = 100;
//		@Option(name="-outputReads",usage=" Outputs one line per read (default false)")
//		protected boolean outputReads = false;
		@Option(name="-useCpgsToFilter",usage=" Use CpGs and CpHs to filter if true, otherwise just CpHs (default false)")
		protected boolean useCpgsToFilter = false;
		@Option(name="-outputHcphs",usage=" Output HCpH cytosines (can't use more than 1 input file)")
		protected boolean outputHcphs = false;
		@Option(name="-minMapQ",usage="minimum mapping quality (default 0)")
		protected int minMapQ = 0;
//		@Option(name="-minHcphCoverage",usage="the minimum number of total reads to include a HCph (default 10)")
//		protected int minHcphCoverage = 10;
//		@Option(name="-minHcphFrac",usage="minimum methylation fraction to include a HCph")
//		protected double minHcphFrac = 0;
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
			new SamToMethyldbOfflineAllCytocineASM().doMain(args);
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
			// Cph counters
			int numCphConvertedWithFilt = 0;
			int numCphTotalWithFilt = 0;
			int numCphConvertedNoFilt = 0;
			int numCphTotalNoFilt = 0;

			// Iterate through chroms
			if (chrs.size()==0) chrs = MethylDbUtils.CHROMS;
			for (final String chr : chrs)
			{
				SortedMap<String,Cytosine> cytosines = new TreeMap<String,Cytosine>();
				PrintWriter outWriter = Cytosine.outputChromToFile(cytosines, prefix, sampleName, chr, asmFlag);

				
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
					int lastBaseSeen = 0;
					int lastPurge = 0;
					
					TreeSet<Integer> allelePosition = allelePos( inputSam, canPurge, chr );
					
					CloseableIterator<SAMRecord> chrIt = inputSam.query(chr, 0, 0, false);
					
					record: while (chrIt.hasNext())
					{
						
						if (canPurge && (lastBaseSeen > (lastPurge+PURGE_INTERVAL)))
						{
							
							//System.err.printf("On base %d, purging everything before %d\n", lastBaseSeen,lastPurge);
							Cytosine.outputCytocinesToFile(outWriter, cytosines.headMap(Integer.toString(lastPurge)), prefix, sampleName, chr, asmFlag);
							
							// Weird, if i just set cytosines to be the tailMap (as in 1, below) garbage collection doesn't actually clean up
							// the old part (just backed by the original data structure). So i actually have to copy it to a new map. (as in 2)
							//(1) cytosines = cytosines.tailMap(new Integer(lastPurge));
							cytosines = new TreeMap<String,Cytosine>(cytosines.tailMap(Integer.toString(lastPurge))); // (2)
							
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
							int seqLen = Math.min(seq.length(), ref.length());
							
							//TreeSet<Integer> AlleleCoord = alleleOfReads(seq, ref, seqLen, onRefCoord, negStrand);
							//Iterator<Integer> AlleleCoordIt = AlleleCoord.iterator();
							
							//if( !AlleleCoord.isEmpty() )
							//	closeAlleleCoord = AlleleCoord.first();
							boolean seqFlag = true;
							int readsStart = samRecord.getUnclippedEnd() <= alignmentS ? samRecord.getUnclippedEnd() : alignmentS;
							int readsEnd = samRecord.getUnclippedEnd() <= alignmentS ? alignmentS : samRecord.getUnclippedEnd();
							if (!allelePosition.subSet(readsStart, readsEnd).isEmpty()){									
								seqFlag = false;
							}
							char A_BaseUpperCase = '0';
							char B_BaseUpperCase = '0';
							
							for (int i = 0; i < seqLen; i++)
							{
								int alleleChromPos = -1;
								char refi = ref.charAt(i);
								char seqi = seq.charAt(i);
								char preBaseRef = PicardUtils.preBaseRef(i, ref);
								char preBaseSeq = PicardUtils.preBaseSeq(i, seq);
								char nextBaseRef = PicardUtils.nextBaseRef(i, ref);
								char nextBaseSeq = PicardUtils.nextBaseSeq(i, seq);
								//System.err.printf("%d\t%d\t%d\t%c\t%c\n",readsStart,readsEnd,i,refi,seqi);
								//if ((i < (seqLen-1)) && PicardUtils.isCytosine(i,ref)) // The last one is too tricky to deal with since we don't know context
								if (PicardUtils.isCytosine(i,ref,false) && PicardUtils.isCytosine(i,seq,true)) // The last one is too tricky to deal with since we don't know context
								{
									boolean isgch = PicardUtils.isGch(i,ref);
									boolean isgcg = PicardUtils.isGcg(i,ref);
									boolean ishcg = PicardUtils.isHcg(i,ref);
									boolean ishch = PicardUtils.isHch(i,ref);
									//boolean iscpg = PicardUtils.isCpg(i,ref);
									
									boolean conv = PicardUtils.isConverted(i,ref,seq);
									
//									if (iscpg && (nextBaseRef != 'G'))
//									{
//										System.err.printf("Cpg status and next base don't match:\nref=%s\nseq=%s\nref_%d=%c, seq_%d=%c\n",
//												ref, seq, i, refi, i, seqi);
//									}
									

									//if (conv && this.useCpgsToFilter || !iscpg) numConverted++;
									if (conv && (this.useCpgsToFilter || (!ishcg & !isgcg))) numConverted++;
									// If this is the first legal one , note it
									if ((convStart==Integer.MAX_VALUE) && (numConverted>=this.minConv) )
									{
										convStart = i;
									}
									
									if(!outputHcphs & ishch){
										
									}
									else{
										if (seqFlag){
											Cytosine cytosine = findOrCreateCytosine(cytosines, onRefCoord, negStrand, preBaseRef, nextBaseRef, A_BaseUpperCase, B_BaseUpperCase, alleleChromPos);
											if (cytosine.getNextBaseRef() == '0') cytosine.setNextBaseRef(nextBaseRef);
											if (cytosine.getPreBaseRef() == '0') cytosine.setPreBaseRef(preBaseRef);
											this.incrementCytosine(cytosine, seqi, i<convStart, preBaseSeq, nextBaseSeq, seqFlag);
										}
										else{
											if (!allelePosition.subSet(readsStart, readsEnd).isEmpty()){									
												Iterator<Integer> allelePosIt = allelePosition.subSet(readsStart, readsEnd).iterator();
												 while (allelePosIt.hasNext()){
													 alleleChromPos = allelePosIt.next();
													 int pos = negStrand ? Math.abs(alleleChromPos - readsEnd) : Math.abs(alleleChromPos - readsStart);
													 A_BaseUpperCase = negStrand ? PicardUtils.revNucleotide(ref.charAt(pos)) : ref.charAt(pos);
													 boolean seqFlag2 = true;
													 if( (ref.charAt(pos) == 'G' && seq.charAt(pos) == 'A') || (ref.charAt(pos) == 'A' && seq.charAt(pos) == 'G')){
														 B_BaseUpperCase = negStrand ? PicardUtils.revNucleotide(ref.charAt(pos)) : seq.charAt(pos); 
														 seqFlag2 = false;
													 }
													 Cytosine cytosine = findOrCreateCytosine(cytosines, onRefCoord, negStrand, preBaseRef, nextBaseRef, A_BaseUpperCase, B_BaseUpperCase, alleleChromPos);
														if (cytosine.getNextBaseRef() == '0') cytosine.setNextBaseRef(nextBaseRef);
														if (cytosine.getPreBaseRef() == '0') cytosine.setPreBaseRef(preBaseRef);
														if (cytosine.getA_BaseUpperCase() == '0') cytosine.setA_BaseUpperCase(A_BaseUpperCase);
														if (cytosine.getB_BaseUpperCase() == '0') cytosine.setB_BaseUpperCase(B_BaseUpperCase);
														this.incrementCytosine(cytosine, seqi, i<convStart, preBaseSeq, nextBaseSeq, seqFlag2
																);
												 }
												 seqFlag = false;
											}
											/*
											int alleleUpperPos = Integer.MAX_VALUE;
											int alleleLowerPos = Integer.MIN_VALUE;
											if(!allelePosition.tailSet(onRefCoord).isEmpty())
												alleleUpperPos = allelePosition.tailSet(onRefCoord).first();
											if(!allelePosition.headSet(onRefCoord).isEmpty())
												alleleLowerPos = allelePosition.headSet(onRefCoord).last();
											if( allelePosition.subSet(readsStart, readsEnd).contains(alleleUpperPos) && allelePosition.subSet(readsStart, readsEnd).contains(alleleLowerPos) ){
												alleleChromPos = Math.abs(alleleUpperPos - onRefCoord) > Math.abs(alleleLowerPos - onRefCoord) ? alleleLowerPos : alleleUpperPos;
											}
											else if (allelePosition.subSet(readsStart, readsEnd).contains(alleleUpperPos) && !allelePosition.subSet(readsStart, readsEnd).contains(alleleLowerPos)){
												alleleChromPos = alleleUpperPos;
											}
											else if (allelePosition.subSet(readsStart, readsEnd).contains(alleleLowerPos) && !allelePosition.subSet(readsStart, readsEnd).contains(alleleUpperPos)){
												alleleChromPos = alleleLowerPos;
											}
											else{
												System.err.printf("error here!");
											}
											int alleleRelativePos = negStrand ? Math.abs(alleleChromPos - readsEnd) : Math.abs(alleleChromPos - readsStart);
											A_BaseUpperCase = ref.charAt(alleleRelativePos);
											boolean converted = PicardUtils.isConverted(alleleRelativePos, ref, seq);
											if (seq.charAt(alleleRelativePos) == A_BaseUpperCase || converted)
												seqFlag = true;
											else
												B_BaseUpperCase = negStrand ? PicardUtils.revNucleotide(seq.charAt(alleleRelativePos)) : seq.charAt(alleleRelativePos);
											//if( A_BaseUpperCase != '0' && (A_BaseUpperCase == 'T' || A_BaseUpperCase == 'C') ){
												System.err.println(ref);
												System.err.println(seq);
												System.err.println(negStrand);
												System.err.println(onRefCoord);
												System.err.println(A_BaseUpperCase);
												System.err.println( alleleChromPos);
												System.err.println(readsStart);
												System.err.println(readsEnd);
												System.err.println(i);
												System.err.println(ref.charAt(alleleRelativePos));
												System.err.println(allelePosition.subSet(readsStart, readsEnd));
												//System.out.printf("%c\t%d\t%c\t%d\t%d\t%d\n", negStrand, onRefCoord, A_BaseUpperCase, alleleChromPos, readsStart, readsEnd);
												System.exit(1);
											}//
											A_BaseUpperCase = negStrand ? PicardUtils.revNucleotide(ref.charAt(alleleRelativePos)) : ref.charAt(alleleRelativePos);
											Cytosine cytosine = findOrCreateCytosine(cytosines, onRefCoord, negStrand, preBaseRef, nextBaseRef, A_BaseUpperCase, B_BaseUpperCase, alleleChromPos);
											if (cytosine.getNextBaseRef() == '0') cytosine.setNextBaseRef(nextBaseRef);
											if (cytosine.getPreBaseRef() == '0') cytosine.setPreBaseRef(preBaseRef);
											if (cytosine.getA_BaseUpperCase() == '0') cytosine.setA_BaseUpperCase(A_BaseUpperCase);
											if (cytosine.getB_BaseUpperCase() == '0') cytosine.setB_BaseUpperCase(B_BaseUpperCase);
											this.incrementCytosine(cytosine, seqi, i<convStart, preBaseSeq, nextBaseSeq, seqFlag);
											seqFlag = false;
											*/
										}
										
									}
									
									

									//if(iscpg)
									if (isgcg || ishcg)
									{
										if (i<convStart)
										{
											// In the non-conversion filter zone
											filteredOutCounter++;
											//System.err.printf("Rec %d\tpos=%d\n",recCounter,i);
										}
										else
										{
											// Past the non-conversion filter, use it
											usedCounter++;
										}

									}
									
									if (isgch || ishch){
										
										numCphTotalNoFilt++;
										if (conv) numCphConvertedNoFilt++;
										
										if (i>=convStart)
										{
											numCphTotalWithFilt++;
											if (conv) numCphConvertedWithFilt++;
										}
									}
									



								} // IsCytosine
/*
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
									Cytosine oppositeCytosine = findOrCreateCytosine(cytosines, onRefCoord, !negStrand, preBaseRefRev, nextBaseRefRev, A_BaseUpperCase, B_BaseUpperCase, alleleChromPos);
							//		
									this.incrementOppositeCytosine(oppositeCytosine, seqi);
								}
*/								


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
				if (cytosines != null) 	Cytosine.outputCytocinesToFile(outWriter, cytosines, prefix, sampleName, chr, asmFlag);
				outWriter.close();

			}
			double frac = (double)filteredOutCounter/((double)usedCounter+(double)filteredOutCounter);
			System.err.printf("Lost %f%% due to non-converion filter\n%d CpGs filtered for non-conversion, %d CpGs used (MinConv=%d,UseCpgs=%s)\n",
					frac*100.0, filteredOutCounter, usedCounter, this.minConv, String.valueOf(this.useCpgsToFilter));
			System.err.printf("Found %d reads total\n", recCounter);
			System.err.printf("CpH conversion rate: before filter=%f, after filter=%f\n",
					100.0 * ((double)numCphConvertedNoFilt/(double)numCphTotalNoFilt),
					100.0 * ((double)numCphConvertedWithFilt/(double)numCphTotalWithFilt));
		}


		protected static Cytosine findOrCreateCytosine(Map<String,Cytosine> cytosines, int onRefCoord, boolean negStrand, char preBaseRef, char nextBaseRef, char A_BaseUpperCase, char B_BaseUpperCase, int alleleChromPos)
		{
			
			String index = Integer.toString(onRefCoord).concat(Integer.toString(alleleChromPos));
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

		
		protected TreeSet<Integer> allelePos( SAMFileReader inputSam, boolean canPurge, String chr )
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
					int readsStart = samRecord.getUnclippedEnd() <= alignmentS ? samRecord.getUnclippedEnd() : alignmentS;
					//int readsEnd = samRecord.getUnclippedEnd() <= alignmentS ? alignmentS : samRecord.getUnclippedEnd();

					if ((this.outputHcphs) && (alignmentS < lastBaseSeen))
					{
						System.err.printf("BAM must be ordered in order to handle -outputCphs: %d<%d\n",alignmentS, lastBaseSeen);
						System.exit(1);
					}
					lastBaseSeen = alignmentS;

					int seqLen = Math.min(seq.length(), ref.length());
					
					for (int i = 0; i < seqLen; i++){
						char refi = ref.charAt(i);
						//char seqi = seq.charAt(i);
						char preBaseRef = PicardUtils.preBaseRef(i, ref);
						char nextBaseRef = PicardUtils.nextBaseRef(i, ref);

						if(negStrand){
							
						}
						else{
							if ((PicardUtils.isGuanine(i,ref) && PicardUtils.isAdenine(i,seq)) || (PicardUtils.isAdenine(i,ref) && PicardUtils.isGuanine(i,seq)) && (nextBaseRef != 'C' && preBaseRef != 'C')){
								Integer readPos = readsStart + i;
								allelePosition.add(readPos);
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
		


}
