package edu.usc.epigenome.scripts;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.SortedMap;
import java.util.TreeMap;

import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.util.CloseableIterator;

//import org.biojava.bio.symbol.Location;
import org.biojava.bio.symbol.LocationTools;
import org.biojava3.genome.parsers.gff.FeatureList;
import org.biojava3.genome.parsers.gff.GFF3Reader;
import org.biojava3.genome.parsers.gff.Location;
import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;

import edu.usc.epigenome.genomeLibs.PicardUtils;

public class SamToSnpBisulfitePrior {
	//final private static String prefix = "methylCGsRich_ASM_AllSnp_";
	final private static String USAGE = "SamToBaseQ [opts] sampleName snpFileName cgiGffFileName bamFilePath ...";

	final private static int PURGE_INTERVAL = 20000; // We purge our stored Cpgs once we get this many bases past them.
	//final private static int ALLELE_GA_NUMBER = 1; //for each reads, when there are more than 1 GA position is different from reference sequence, we define it belongs to another allele. 
	
	/**
	 * @param args
	 */
	
	@Option(name="-minMapQ",usage="minimum mapping quality (default 30)")
	protected int minMapQ = 30;
	@Option(name="-minReadCov",usage="minimum read coverage (default 5)")
	protected int minReadCov = 5;
	@Option(name="-minBaseQual",usage="minimum Base quality (default 30)")
	protected int minBaseQual = 30;
	@Option(name="-snpType",usage="snp position type: 1: in cytosine; 2: next base of cytosine; 3: 3rd base of cytosine (default: 1)")
	protected int snpType = 1;
	@Option(name="-debug",usage=" Debugging statements (default false)")
	protected boolean debug = false;
	// receives other command line parameters than options
	@Argument
	private List<String> stringArgs = new ArrayList<String>();
	
	private class SnpData{
		private String snpName = null;
		private String snpChr = null;
		private long snpPos = 0;
		private String snpAlleleA = null;
		private String snpAlleleB = null;
		private String snpAlleleAstr = null;
		private String snpAlleleBstr = null;
		private String snpAlleleAstr3Char = null;
		private String snpAlleleBstr3Char = null;
		private int numC = 0;
		private int numT = 0;
		private int numTotal = 0;
		private boolean cgiStat = false;
		//private double perC = 0;
		//private double perT = 0;
		
	}
	
	
	
	public static void main(String[] args) throws Exception
	{
		new SamToSnpBisulfitePrior().doMain(args);
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
		String cgiFileName = stringArgs.remove(0);
		String bamFilePath = stringArgs.remove(0);
		BufferedReader br = new BufferedReader(new FileReader(snpFileName));
		//BufferedReader cgiBr = new BufferedReader(new FileReader(cgiFileName));
		//SortedMap<String,Location> cgi = new TreeMap<String,Location>();
		FeatureList cgi = GFF3Reader.read(cgiFileName);
		//cgiLocation(cgiFileName, cgi);
		String line;
		
		String preChr = "chr1";
		String fn = sampleName + "Merge_chr1.NODUPS.sorted.calmd.NODUPS.bam";
		File inputSamOrBamFile = new File(bamFilePath,fn);
		SAMFileReader inputSam = new SAMFileReader(inputSamOrBamFile);
		inputSam.setValidationStringency(SAMFileReader.ValidationStringency.SILENT);
		String fnSnp = sampleName + "_SNP." + minMapQ + "-" + minBaseQual + "-" + minReadCov + "-" +  "SNP"+ snpType + "." + "txt";
		//String fnSnp2 = sampleName + "_SNP." + minMapQ + "-" + minBaseQual + "-" + minReadCov + "-" +  "SNP2"  + "." + "txt";
		//String fnSnp3 = sampleName + "_SNP." + minMapQ + "-" + minBaseQual + "-" + minReadCov + "-" +  "SNP3"  + "." + "txt";
		
		PrintWriter snpWriter = new PrintWriter(new File(fnSnp));
		//PrintWriter snpWriter2 = new PrintWriter(new File(fnSnp2));
		//PrintWriter snpWriter3 = new PrintWriter(new File(fnSnp3));
		
		//using nearest cytosine's location to be the key, for type 1, there is only one key for each snp, but for type2 & 3, it is possible to have 2 cytosine in both direction. 
		TreeMap<Integer,SnpData> snpSet = new TreeMap<Integer,SnpData>();
		
		while( (line = br.readLine()) != null){
			String[] tmpArray = line.split("\t");
			String chr = "chr" + tmpArray[1];
			//boolean xFlag = false;
			if(chr.equalsIgnoreCase("chrX") || chr.equalsIgnoreCase("chrY") || chr.equalsIgnoreCase("chrM")){
			//	xFlag = true;
				continue;
			}
			if(!chr.equalsIgnoreCase(preChr)){
				inputSam.close();
				snpSetWriter(snpSet,snpWriter);
				snpSet.clear();
				fn = sampleName + "Merge_" + chr + ".NODUPS.sorted.calmd.NODUPS.bam";
				File inputNewSamOrBamFile = new File(bamFilePath,fn);
				System.err.println(inputNewSamOrBamFile.getName());
				inputSam = new SAMFileReader(inputNewSamOrBamFile);
				inputSam.setValidationStringency(SAMFileReader.ValidationStringency.SILENT);
				preChr = chr;
			}
			
			Integer snpPosition =(int) Math.round(Double.parseDouble(tmpArray[2]));
			Integer cytosineLocation = null;
			Integer cytosineLocationBack = null;
			Location snpLocation = new Location(snpPosition, snpPosition);
			SnpData snpOutput = new SnpData();
			SnpData snpOutputBack = new SnpData();
			snpOutput.snpChr = chr;
			snpOutput.snpPos = snpPosition;
			if(!cgi.selectOverlapping(chr, snpLocation,true).isEmpty()){
				snpOutput.cgiStat = true;
			}
			snpOutput.snpName = tmpArray[0];
			snpOutput.snpAlleleA = tmpArray[3];
			snpOutput.snpAlleleB = tmpArray[4];
			snpOutput.snpAlleleAstr = tmpArray[5];
			snpOutput.snpAlleleBstr = tmpArray[6];
			
			//boolean effectReadInNegStrand = false;
			
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
					byte refi = (byte) ref.charAt(i);
					if(snpType == 1){
						if( !tmpArray[5].equalsIgnoreCase("C") && !tmpArray[5].equalsIgnoreCase("G") && !tmpArray[6].equalsIgnoreCase("C") && !tmpArray[6].equalsIgnoreCase("G") ){
							break record;
						}
						cytosineLocation = snpPosition;
						//if((PicardUtils.isCytosine(i, ref, false) && negStrand) || PicardUtils.revNucleotide(ref.charAt(i)){
							//effectReadInNegStrand = true;
						//}
						//else{
							
					//	}
						/*
							if( tmpArray[5].equalsIgnoreCase(tmpArray[6]) ){
								if(ref.charAt(i) == 'C'){
										if(i < seqLen-2){
											if(snpOutput.snpRefstr3Char.isEmpty()){
												snpOutput.snpRefstr3Char = ref.substring(i,i+2);
												if(ref.charAt(i) != seq.charAt(i) && !PicardUtils.isCytosine(i, seq, true)){
													snpOutput.snpAllelestr3Char = seq.substring(i,i+2);
												}
												else{
													snpOutput.snpAllelestr3Char = snpOutput.snpRefstr3Char;
												}
											}
											else if(ref.charAt(i) != seq.charAt(i) && !PicardUtils.isCytosine(i, seq, true) && snpOutput.snpRefstr3Char.equalsIgnoreCase(snpOutput.snpAllelestr3Char)){
												snpOutput.snpAllelestr3Char = seq.substring(i,i+2);
											}
										}
								}
								else if(ref.charAt(i) == 'G'){
									if(i > 1){
										if(snpOutput.snpRefstr3Char.isEmpty()){
											snpOutput.snpRefstr3Char = "" + PicardUtils.revNucleotide(ref.charAt(i)) + PicardUtils.revNucleotide(ref.charAt(i-1)) + PicardUtils.revNucleotide(ref.charAt(i-2));
											snpOutput.snpAllelestr3Char = snpOutput.snpRefstr3Char;
										}
										else if(ref.charAt(i) != seq.charAt(i) && snpOutput.snpRefstr3Char.equalsIgnoreCase(snpOutput.snpAllelestr3Char)){
											snpOutput.snpAllelestr3Char = "" + PicardUtils.revNucleotide(seq.charAt(i)) + PicardUtils.revNucleotide(seq.charAt(i-1)) + PicardUtils.revNucleotide(seq.charAt(i-2));
											
										}
									}
								}
							}
							else{
								if(ref.charAt(i) == 'C'){
									if(i < seqLen-2){
										if(snpOutput.snpRefstr3Char.isEmpty()){
											snpOutput.snpRefstr3Char = ref.substring(i,i+2);
											if(ref.charAt(i) != seq.charAt(i) && !PicardUtils.isCytosine(i, seq, true)){
												snpOutput.snpAllelestr3Char = seq.substring(i,i+2);
											}
											else{
												snpOutput.snpAllelestr3Char = snpOutput.snpRefstr3Char;
											}
										}
										else if(ref.charAt(i) != seq.charAt(i) && !PicardUtils.isCytosine(i, seq, true) && snpOutput.snpRefstr3Char.equalsIgnoreCase(snpOutput.snpAllelestr3Char)){
											snpOutput.snpAllelestr3Char = seq.substring(i,i+2);
										}
									}
								}
								else if(ref.charAt(i) != 'C' && (ref.charAt(i) == tmpArray[5].charAt(0) || ref.charAt(i) == tmpArray[6].charAt(0))){
									if(snpOutput.snpRefstr3Char.isEmpty()){
										snpOutput.snpRefstr3Char = ref.substring(i,i+2);
										if(ref.charAt(i) != seq.charAt(i)){
											snpOutput.snpAllelestr3Char = seq.substring(i,i+2);
										}
										else{
											snpOutput.snpAllelestr3Char = snpOutput.snpRefstr3Char;
										}
									}
									else if(ref.charAt(i) != seq.charAt(i) && snpOutput.snpRefstr3Char.equalsIgnoreCase(snpOutput.snpAllelestr3Char)){
										snpOutput.snpAllelestr3Char = seq.substring(i,i+2);
									}
								}
								else if(ref.charAt(i) != 'C' && (ref.charAt(i) != tmpArray[5].charAt(0) && ref.charAt(i) != tmpArray[6].charAt(0))){
									if(snpOutput.snpRefstr3Char.isEmpty()){
										snpOutput.snpRefstr3Char = ref.substring(i,i+2);
										if(ref.charAt(i) != seq.charAt(i)){
											snpOutput.snpAllelestr3Char = seq.substring(i,i+2);
										}
										else{
											snpOutput.snpAllelestr3Char = snpOutput.snpRefstr3Char;
										}
									}
									else if(ref.charAt(i) != seq.charAt(i) && snpOutput.snpRefstr3Char.equalsIgnoreCase(snpOutput.snpAllelestr3Char)){
										snpOutput.snpAllelestr3Char = seq.substring(i,i+2);
									}
								}
							}
							*/
							if(ref.charAt(i) == 'C'){
								if(i < seqLen-2){
									if(snpOutput.snpAlleleAstr3Char.isEmpty() || snpOutput.snpAlleleBstr3Char.isEmpty()){
										if(PicardUtils.revNucleotide(tmpArray[5].charAt(0)) == 'C' || PicardUtils.revNucleotide(tmpArray[6].charAt(0)) == 'C'){
											snpOutput.snpAlleleAstr3Char = PicardUtils.revNucleotide(tmpArray[5].charAt(0)) + ref.substring(i+1,i+2);
											snpOutput.snpAlleleBstr3Char = PicardUtils.revNucleotide(tmpArray[6].charAt(0)) + ref.substring(i+1,i+2);
										}
										else{
											snpOutput.snpAlleleAstr3Char = tmpArray[5] + ref.substring(i+1,i+2);
											snpOutput.snpAlleleBstr3Char = tmpArray[6] + ref.substring(i+1,i+2);
										}
										
									}
								}
								snpOutput.numTotal++;
							}
							else if(PicardUtils.revNucleotide(ref.charAt(i)) == 'C'){
								if(i > 1){
									if(snpOutput.snpAlleleAstr3Char.isEmpty() || snpOutput.snpAlleleBstr3Char.isEmpty()){
										if(PicardUtils.revNucleotide(tmpArray[5].charAt(0)) == 'C' || PicardUtils.revNucleotide(tmpArray[6].charAt(0)) == 'C'){
											snpOutput.snpAlleleAstr3Char = "" + PicardUtils.revNucleotide(tmpArray[5].charAt(0)) + PicardUtils.revNucleotide(ref.charAt(i-1)) + PicardUtils.revNucleotide(ref.charAt(i-2));
											snpOutput.snpAlleleBstr3Char = "" + PicardUtils.revNucleotide(tmpArray[6].charAt(0)) + PicardUtils.revNucleotide(ref.charAt(i-1)) + PicardUtils.revNucleotide(ref.charAt(i-2));
										}
										else{
											snpOutput.snpAlleleAstr3Char = tmpArray[5] + PicardUtils.revNucleotide(ref.charAt(i-1)) + PicardUtils.revNucleotide(ref.charAt(i-2));
											snpOutput.snpAlleleBstr3Char = tmpArray[6] + PicardUtils.revNucleotide(ref.charAt(i-1)) + PicardUtils.revNucleotide(ref.charAt(i-2));
										}
										
									}
								}
								if((tmpArray[5].charAt(0) == 'C' && tmpArray[6].charAt(0) == 'G') || (tmpArray[6].charAt(0) == 'C' && tmpArray[5].charAt(0) == 'G')){
									snpOutput.numTotal++;
								}
							}
							else{
								if(ref.charAt(i) == tmpArray[5].charAt(0) || ref.charAt(i) == tmpArray[6].charAt(0)){
									if(PicardUtils.revNucleotide(tmpArray[5].charAt(0)) == 'C' || PicardUtils.revNucleotide(tmpArray[6].charAt(0)) == 'C'){
										if(i > 1){
											if(snpOutput.snpAlleleAstr3Char.isEmpty() || snpOutput.snpAlleleBstr3Char.isEmpty()){
												snpOutput.snpAlleleAstr3Char = "" + PicardUtils.revNucleotide(tmpArray[5].charAt(0)) + PicardUtils.revNucleotide(ref.charAt(i-1)) + PicardUtils.revNucleotide(ref.charAt(i-2));
												snpOutput.snpAlleleBstr3Char = "" + PicardUtils.revNucleotide(tmpArray[6].charAt(0)) + PicardUtils.revNucleotide(ref.charAt(i-1)) + PicardUtils.revNucleotide(ref.charAt(i-2));
											}
										}
									}
									else{
										if(i < seqLen-2){
											if(snpOutput.snpAlleleAstr3Char.isEmpty() || snpOutput.snpAlleleBstr3Char.isEmpty()){
												snpOutput.snpAlleleAstr3Char = tmpArray[5] + ref.substring(i+1,i+2);
												snpOutput.snpAlleleBstr3Char = tmpArray[6] + ref.substring(i+1,i+2);
											}
										}
										snpOutput.numTotal++;
									}
								}
								else{
									if(tmpArray[5].charAt(0) == 'C' || tmpArray[6].charAt(0) == 'C'){
										if(i > 1){
											if(snpOutput.snpAlleleAstr3Char.isEmpty() || snpOutput.snpAlleleBstr3Char.isEmpty()){
												snpOutput.snpAlleleAstr3Char = tmpArray[5] + PicardUtils.revNucleotide(ref.charAt(i-1)) + PicardUtils.revNucleotide(ref.charAt(i-2));
												snpOutput.snpAlleleBstr3Char = tmpArray[6] + PicardUtils.revNucleotide(ref.charAt(i-1)) + PicardUtils.revNucleotide(ref.charAt(i-2));
											}
										}
									}
									else{
										if(i < seqLen-2){
											if(snpOutput.snpAlleleAstr3Char.isEmpty() || snpOutput.snpAlleleBstr3Char.isEmpty()){
												snpOutput.snpAlleleAstr3Char = PicardUtils.revNucleotide(tmpArray[5].charAt(0)) + ref.substring(i+1,i+2);
												snpOutput.snpAlleleBstr3Char = PicardUtils.revNucleotide(tmpArray[6].charAt(0)) + ref.substring(i+1,i+2);
											}
										}
										snpOutput.numTotal++;
									}
								}
							}

						if(PicardUtils.isCytosine(i, seq, false)){
							snpOutput.numC++;
						}
						else if(PicardUtils.isThymine(i, seq)){
							snpOutput.numT++;
						}
						 
						
					}
					else if(snpType == 2){		
						if(i >= 1 && i < seqLen-2){
								// two direction both have cytosine
								if(PicardUtils.isCytosine(i-1, ref, false) && PicardUtils.revNucleotide(ref.charAt(i+1)) == 'C'){
									cytosineLocation = snpPosition-1;
									cytosineLocationBack = snpPosition+1;
									if(snpOutputBack.snpAlleleA.isEmpty()){
										snpOutputBack.snpChr = chr;
										snpOutputBack.snpPos = snpPosition;
										if(!cgi.selectOverlapping(chr, snpLocation,true).isEmpty()){
											snpOutputBack.cgiStat = true;
										}
										snpOutputBack.snpName = tmpArray[0];
										snpOutputBack.snpAlleleA = tmpArray[3];
										snpOutputBack.snpAlleleB = tmpArray[4];
										snpOutputBack.snpAlleleAstr = tmpArray[5];
										snpOutputBack.snpAlleleBstr = tmpArray[6];
									}
									if(!negStrand){
										if(snpOutput.snpAlleleAstr3Char.isEmpty() || snpOutput.snpAlleleBstr3Char.isEmpty()){
											if( ref.charAt(i) == tmpArray[5].charAt(0) || ref.charAt(i) == tmpArray[6].charAt(0) ){
												snpOutput.snpAlleleAstr3Char = "" + ref.charAt(i-1) + tmpArray[5].charAt(0) + ref.charAt(i+1);
												snpOutput.snpAlleleBstr3Char =  "" + ref.charAt(i-1) + tmpArray[6].charAt(0) + ref.charAt(i+1);
											}
											else{
												snpOutput.snpAlleleAstr3Char = "" + ref.charAt(i-1) + PicardUtils.revNucleotide(tmpArray[5].charAt(0)) + ref.charAt(i+1);
												snpOutput.snpAlleleBstr3Char =  "" + ref.charAt(i-1) + PicardUtils.revNucleotide(tmpArray[6].charAt(0)) + ref.charAt(i+1);
											}		
										}
										if(snpOutputBack.snpAlleleAstr3Char.isEmpty() || snpOutputBack.snpAlleleBstr3Char.isEmpty()){
											if( ref.charAt(i) == tmpArray[5].charAt(0) || ref.charAt(i) == tmpArray[6].charAt(0) ){
												snpOutputBack.snpAlleleAstr3Char = "" + PicardUtils.revNucleotide(ref.charAt(i+1)) + PicardUtils.revNucleotide(tmpArray[5].charAt(0)) + PicardUtils.revNucleotide(ref.charAt(i-1));
												snpOutputBack.snpAlleleBstr3Char =  "" + PicardUtils.revNucleotide(ref.charAt(i+1)) + PicardUtils.revNucleotide(tmpArray[6].charAt(0)) + PicardUtils.revNucleotide(ref.charAt(i-1));
											}
											else{
												snpOutputBack.snpAlleleAstr3Char = "" + PicardUtils.revNucleotide(ref.charAt(i+1)) + tmpArray[5].charAt(0) + PicardUtils.revNucleotide(ref.charAt(i-1));
												snpOutputBack.snpAlleleBstr3Char =  "" + PicardUtils.revNucleotide(ref.charAt(i+1)) + tmpArray[6].charAt(0) + PicardUtils.revNucleotide(ref.charAt(i-1));
											}
										}
										if(PicardUtils.isCytosine(i-1, seq, false)){
											snpOutput.numC++;
										}
										else if(PicardUtils.isThymine(i-1, seq)){
											snpOutput.numT++;
										}
												
									}
									else{
										if(snpOutput.snpAlleleAstr3Char.isEmpty() || snpOutput.snpAlleleBstr3Char.isEmpty()){
											if( ref.charAt(i) == tmpArray[5].charAt(0) || ref.charAt(i) == tmpArray[6].charAt(0) ){
												snpOutput.snpAlleleAstr3Char = "" + PicardUtils.revNucleotide(ref.charAt(i+1)) + PicardUtils.revNucleotide(tmpArray[5].charAt(0)) + PicardUtils.revNucleotide(ref.charAt(i-1));
												snpOutput.snpAlleleBstr3Char =  "" + PicardUtils.revNucleotide(ref.charAt(i+1)) + PicardUtils.revNucleotide(tmpArray[6].charAt(0)) + PicardUtils.revNucleotide(ref.charAt(i-1));
											}
											else{
												snpOutput.snpAlleleAstr3Char = "" + PicardUtils.revNucleotide(ref.charAt(i+1)) + tmpArray[5].charAt(0) + PicardUtils.revNucleotide(ref.charAt(i-1));
												snpOutput.snpAlleleBstr3Char =  "" + PicardUtils.revNucleotide(ref.charAt(i+1)) + tmpArray[6].charAt(0) + PicardUtils.revNucleotide(ref.charAt(i-1));
											}		
										}
										if(snpOutputBack.snpAlleleAstr3Char.isEmpty() || snpOutputBack.snpAlleleBstr3Char.isEmpty()){
											if( ref.charAt(i) == tmpArray[5].charAt(0) || ref.charAt(i) == tmpArray[6].charAt(0) ){
												snpOutputBack.snpAlleleAstr3Char = "" + ref.charAt(i-1) + tmpArray[5].charAt(0) + ref.charAt(i+1);
												snpOutputBack.snpAlleleBstr3Char =  "" + ref.charAt(i-1) + tmpArray[6].charAt(0) + ref.charAt(i+1);
											}
											else{
												snpOutputBack.snpAlleleAstr3Char = "" + ref.charAt(i-1) + PicardUtils.revNucleotide(tmpArray[5].charAt(0)) + ref.charAt(i+1);
												snpOutputBack.snpAlleleBstr3Char =  "" + ref.charAt(i-1) + PicardUtils.revNucleotide(tmpArray[6].charAt(0)) + ref.charAt(i+1);
											}
										}
										if(PicardUtils.isCytosine(i-1, seq, false)){
											snpOutputBack.numC++;
										}
										else if(PicardUtils.isThymine(i-1, seq)){
											snpOutputBack.numT++;
										}
									}
									continue record;
								}
								else{
									
									if(snpOutput.snpAlleleAstr3Char.isEmpty() || snpOutput.snpAlleleBstr3Char.isEmpty()){
										if(PicardUtils.isCytosine(i-1, ref, false) ||  PicardUtils.revNucleotide(ref.charAt(i+1)) == 'C'){
											if(PicardUtils.isCytosine(i-1, ref, false)){
												if( ref.charAt(i) == tmpArray[5].charAt(0) || ref.charAt(i) == tmpArray[6].charAt(0) ){
													snpOutput.snpAlleleAstr3Char = "" + ref.charAt(i-1) + tmpArray[5].charAt(0) + ref.charAt(i+1);
													snpOutput.snpAlleleBstr3Char =  "" + ref.charAt(i-1) + tmpArray[6].charAt(0) + ref.charAt(i+1);
												}
												else{
													snpOutput.snpAlleleAstr3Char = "" + ref.charAt(i-1) + PicardUtils.revNucleotide(tmpArray[5].charAt(0)) + ref.charAt(i+1);
													snpOutput.snpAlleleBstr3Char =  "" + ref.charAt(i-1) + PicardUtils.revNucleotide(tmpArray[6].charAt(0)) + ref.charAt(i+1);
												}
												if(negStrand){
													cytosineLocation = snpPosition+1;
												}
												else{
													cytosineLocation = snpPosition-1;
												}
											}
											else{
												if( ref.charAt(i) == tmpArray[5].charAt(0) || ref.charAt(i) == tmpArray[6].charAt(0) ){
													snpOutput.snpAlleleAstr3Char = "" + PicardUtils.revNucleotide(ref.charAt(i+1)) + PicardUtils.revNucleotide(tmpArray[5].charAt(0)) + PicardUtils.revNucleotide(ref.charAt(i-1));
													snpOutput.snpAlleleBstr3Char =  "" + PicardUtils.revNucleotide(ref.charAt(i+1)) + PicardUtils.revNucleotide(tmpArray[6].charAt(0)) + PicardUtils.revNucleotide(ref.charAt(i-1));
												}
												else{
													snpOutput.snpAlleleAstr3Char = "" + PicardUtils.revNucleotide(ref.charAt(i+1)) + tmpArray[5].charAt(0) + PicardUtils.revNucleotide(ref.charAt(i-1));
													snpOutput.snpAlleleBstr3Char =  "" + PicardUtils.revNucleotide(ref.charAt(i+1)) + tmpArray[6].charAt(0) + PicardUtils.revNucleotide(ref.charAt(i-1));
												}
												if(negStrand){
													cytosineLocation = snpPosition-1;
												}
												else{
													cytosineLocation = snpPosition+1;
												}
											}
											
										}
									}
								}
							}
						
						
						if(i >= 1){
							if(PicardUtils.isCytosine(i-1, seq, false)){
								snpOutput.numC++;
							}
							else if(PicardUtils.isThymine(i-1, seq)){
								snpOutput.numT++;
							}
						}
						
					}
					else if(snpType == 3){
						if(i >= 2 && i < seqLen-3){
								if(snpOutput.snpAlleleAstr3Char.isEmpty() || snpOutput.snpAlleleBstr3Char.isEmpty()){
									if(PicardUtils.isCytosine(i-2, ref, false) ||  PicardUtils.revNucleotide(ref.charAt(i+2)) == 'C'){
										if(PicardUtils.isCytosine(i-2, ref, false)){
											if( ref.charAt(i+1) == tmpArray[5].charAt(0) || ref.charAt(i+1) == tmpArray[6].charAt(0) ){
												snpOutput.snpAlleleAstr3Char = "" + ref.charAt(i-1) + ref.charAt(i) + tmpArray[5].charAt(0);
												snpOutput.snpAlleleBstr3Char =  "" + ref.charAt(i-1) + ref.charAt(i) + tmpArray[6].charAt(0);
											}
											else{
												snpOutput.snpAlleleAstr3Char = "" + ref.charAt(i-1) + ref.charAt(i) + PicardUtils.revNucleotide(tmpArray[5].charAt(0));
												snpOutput.snpAlleleBstr3Char =  "" + ref.charAt(i-1) + ref.charAt(i) + PicardUtils.revNucleotide(tmpArray[6].charAt(0));
											}
											if(negStrand){
												cytosineLocation = snpPosition+2;
											}
											else{
												cytosineLocation = snpPosition-2;
											}
										}
										else{
											if( ref.charAt(i+1) == tmpArray[5].charAt(0) || ref.charAt(i+1) == tmpArray[6].charAt(0) ){
												snpOutput.snpAlleleAstr3Char = "" + PicardUtils.revNucleotide(tmpArray[5].charAt(0)) + PicardUtils.revNucleotide(ref.charAt(i)) + PicardUtils.revNucleotide(ref.charAt(i-1));
												snpOutput.snpAlleleBstr3Char =  "" + PicardUtils.revNucleotide(tmpArray[6].charAt(0)) + PicardUtils.revNucleotide(ref.charAt(i)) + PicardUtils.revNucleotide(ref.charAt(i-1));
											}
											else{
												snpOutput.snpAlleleAstr3Char = "" + tmpArray[5].charAt(0) + PicardUtils.revNucleotide(ref.charAt(i)) + PicardUtils.revNucleotide(ref.charAt(i-1));
												snpOutput.snpAlleleBstr3Char =  "" + tmpArray[6].charAt(0) + PicardUtils.revNucleotide(ref.charAt(i)) + PicardUtils.revNucleotide(ref.charAt(i-1));
											}
											if(negStrand){
												cytosineLocation = snpPosition-2;
											}
											else{
												cytosineLocation = snpPosition+2;
											}
										}
										
									}
								}
						}

						if(i >= 2){
							if(PicardUtils.isCytosine(i-2, seq, false)){
								snpOutput.numC++;
							}
							else if(PicardUtils.isThymine(i-2, seq)){
								snpOutput.numT++;
							}
						}
					
					}
					else{
						System.err.println("error! no snpType definition!");
						System.exit(1);
					}
					
				}
				catch (Exception e)
				{
					System.err.println("-----------------------------------------");
					System.err.println("Couldn't handle seq #");
					System.err.println(seq);
					e.printStackTrace(System.err);
					System.err.println("-----------------------------------------");
//					chrIt.close();
//					System.exit(1);
				}
				
			}// one sam record over
			chrIt.close();
			if(snpType == 2){
				snpSet.put(cytosineLocation, snpOutput);
				if(!snpOutputBack.snpAlleleAstr3Char.isEmpty() && cytosineLocationBack != null){
					snpSet.put(cytosineLocationBack, snpOutputBack);
				}
			}
			else{
				snpSet.put(cytosineLocation, snpOutput);
			}
		}//one line of snp over
	}
	
//	private static void cgiLocation(String cgiFileNameGff, TreeMap<String,Location> cgiLocation) throws IOException{

//		FeatureList cgiFeatures = GFF3Reader.read(cgiFileNameGff);
//		cgi = cgiFeatures.bounds();
//	}

	private static void snpSetWriter(TreeMap<Integer,SnpData> snpSet, PrintWriter snpWriter){
		Iterator<Integer> snpSetIt = snpSet.keySet().iterator();
		while(snpSetIt.hasNext()){
			SnpData snp = snpSet.get(snpSetIt.next());
			snpWriter.println(snp.snpName + "\t" + snp.snpChr + "\t" + snp.snpPos + "\t" + snp.snpAlleleA + "\t" + snp.snpAlleleB + "\t" + snp.snpAlleleAstr + "\t" + snp.snpAlleleBstr + "\t" + snp.snpAlleleAstr3Char + "\t" + snp.snpAlleleBstr3Char + "\t" + snp.numC + "\t" + snp.numT + "\t" + snp.numTotal);
			
		}
	}
}
