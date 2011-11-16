package edu.usc.epigenome.scripts;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.TreeMap;

import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;

public class OverlapAsmBlockCpgReadsCriteria {

	/**
	 * @param args
	 */

	final private static String USAGE = "OverlapAsmBlockCpgReadsCriteria [opt] ASMBlockFileDirecotryPath";
	final private static String prefix = "methylCGsRich_ASM_AllSnp_";
	
	@Option(name="-chrom",multiValued=true,usage="One or more chroms, eg. --chrom chr1 --chrom chr5")
	protected List<String> chrs = new ArrayList<String>(25);
	@Option(name="-minAlleleCpgReads",usage="minimum number in one ASM Block, default value is: 4")
	protected int minAlleleCpgReads = 5;
	@Option(name="-pValue",usage="adjusted p value thresh hold for significant ASM block, default value is: 0.01")
	protected double pValue = 0.01;
	@Option(name="-methyRatio",usage="methylation ratio thresh hold for significant ASM block, default value is: 2")
	protected double methyRatio = 2;
	@Option(name="-alterpValue",usage="unjusted p value thresh hold for not significant ASM block, default value is: 0.05")
	protected double alterpValue = 0.05;
	@Option(name="-alterMethyRatio",usage="methylation ratio thresh hold for not significant ASM block, default value is: 1.5")
	protected double alterMethyRatio = 1.5;
	
	@Argument
	private List<String> stringArgs = new ArrayList<String>();
	
	public static void main(String[] args)
	throws Exception
	{
		new OverlapAsmBlockCpgReadsCriteria().doMain(args);
	}

	public void doMain(String[] args)
	throws Exception {
		CmdLineParser parser = new CmdLineParser(this);
		parser.setUsageWidth(80);
		try
		{
			parser.parseArgument(args);
			if (stringArgs.size() < 1) throw new CmdLineException(USAGE);
			
			
		}
		catch (CmdLineException e)
		{
			System.err.println(e.getMessage());
			parser.printUsage(System.err);
			System.err.println();
			return;
		}


		
		
		String snpFilePath = stringArgs.remove(0);

		

		String normalFn = prefix + "normalMerge_pValue_CpGblock." + minAlleleCpgReads + ".txt";
		String tumorFn = prefix + "tumorMerge_pValue_CpGblock." + minAlleleCpgReads + ".txt";
		String normalGtfFn = prefix + "normalMerge_CpGblock." + minAlleleCpgReads + ".gtf";
		String tumorGtfFn = prefix + "tumorMerge_CpGblock." + minAlleleCpgReads + ".gtf";
		PrintWriter normalWriter = new PrintWriter(new File(normalFn));
		PrintWriter tumorWriter = new PrintWriter(new File(tumorFn));
		PrintWriter normalGtfWriter = new PrintWriter(new File(normalGtfFn));
		PrintWriter tumorGtfWriter = new PrintWriter(new File(tumorGtfFn));
		String normalOverlapFn = prefix + "normalMerge_pValue_CpGblock_overlap." + minAlleleCpgReads + ".txt";
		String tumorOverlapFn = prefix + "tumorMerge_pValue_CpGblock_overlap." + minAlleleCpgReads + ".txt";
		String overlapGtfFn = prefix + "normal-tumorMerge_pValue_CpGblock_overlap." + minAlleleCpgReads + ".gtf";
		PrintWriter normalOverlapWriter = new PrintWriter(new File(normalOverlapFn));
		PrintWriter tumorOverlapWriter = new PrintWriter(new File(tumorOverlapFn));
		PrintWriter overlapGtfWriter = new PrintWriter(new File(overlapGtfFn));
		
		String normalPvalueFn = prefix + "normalMerge_CpGblock.PvalueAdjusted." + minAlleleCpgReads + ".bedGraph";
		String tumorPvalueFn = prefix + "tumorMerge_CpGblock.PvalueAdjusted." + minAlleleCpgReads + ".bedGraph";
		PrintWriter normalPvalueWriter = new PrintWriter(new File(normalPvalueFn));
		PrintWriter tumorPvalueWriter = new PrintWriter(new File(tumorPvalueFn));
		normalPvalueWriter.println("track type=bedGraph name=\"-Log10 adjusted P value normal CpG blocks \" description=\"-Log10 adjusted P value normal CpG blocks\" visibility=full color=200,100,0 altColor=0,100,200 priority=20");
		tumorPvalueWriter.println("track type=bedGraph name=\"-Log10 adjusted P value tumor CpG blocks \" description=\"-Log10 adjusted P value tumor CpG blocks\" visibility=full color=200,100,0 altColor=0,100,200 priority=20");
		String normalPvalueCpGblockOverlapFn = prefix + "normalMerge_CpGblock_overlap.PvalueAdjusted." + minAlleleCpgReads + ".bedGraph";
		String tumorPvalueCpGblockOverlapFn = prefix + "tumorMerge_CpGblock_overlap.PvalueAdjusted." + minAlleleCpgReads + ".bedGraph";
		PrintWriter normalPvalueCpGblockOverlapWriter = new PrintWriter(new File(normalPvalueCpGblockOverlapFn));
		PrintWriter tumorPvalueCpGblockOverlapWriter = new PrintWriter(new File(tumorPvalueCpGblockOverlapFn));
		normalPvalueWriter.println("track type=bedGraph name=\"-Log10 adjusted P value normal CpG blocks (Overlap) \" description=\"-Log10 adjusted P value normal CpG blocks (Overlap)\" visibility=full color=200,100,0 altColor=0,100,200 priority=20");
		tumorPvalueWriter.println("track type=bedGraph name=\"-Log10 adjusted P value tumor CpGblocks (Overlap) \" description=\"-Log10 adjusted P value tumor CpG blocks (Overlap)\" visibility=full color=200,100,0 altColor=0,100,200 priority=20");
		String normalSigPvalueCpGblockOverlapFn = prefix + "normalMerge_sigASM_CpGblock_overlap.PvalueAdjusted." + pValue + ".MethyRatio." + methyRatio + "." + minAlleleCpgReads + ".bedGraph";
		String tumorSigPvalueCpGblockOverlapFn = prefix + "tumorMerge_sigASM_CpGblock_overlap.PvalueAdjusted." + pValue + ".MethyRatio." + methyRatio + "." + minAlleleCpgReads + ".bedGraph";
		PrintWriter normalSigPvalueCpGblockOverlapWriter = new PrintWriter(new File(normalSigPvalueCpGblockOverlapFn));
		PrintWriter tumorSigPvalueCpGblockOverlapWriter = new PrintWriter(new File(tumorSigPvalueCpGblockOverlapFn));
		normalSigPvalueCpGblockOverlapWriter.println("track type=bedGraph name=\"-Log10 adjusted P value normal ASM CpG blocks (Overlap) \" description=\"-Log10 adjusted P value normal ASM CpG blocks (Overlap)\" visibility=full color=200,100,0 altColor=0,100,200 priority=20");
		tumorSigPvalueCpGblockOverlapWriter.println("track type=bedGraph name=\"-Log10 adjusted P value tumor ASM CpGblocks (Overlap) \" description=\"-Log10 adjusted P value tumor ASM CpG blocks (Overlap)\" visibility=full color=200,100,0 altColor=0,100,200 priority=20");
		
		
		String normalPvalueUnjustedFn = prefix + "normalMerge_CpGblock.PvalueUnjusted." + minAlleleCpgReads + ".bedGraph";
		String tumorPvalueUnjustedFn = prefix + "tumorMerge_CpGblock.PvalueUnjusted." + minAlleleCpgReads + ".bedGraph";
		PrintWriter normalPvalueUnjustedWriter = new PrintWriter(new File(normalPvalueUnjustedFn));
		PrintWriter tumorPvalueUnjustedWriter = new PrintWriter(new File(tumorPvalueUnjustedFn));
		normalPvalueUnjustedWriter.println("track type=bedGraph name=\"-Log10 unjusted P value normal CpG blocks \" description=\"-Log10 unjusted P value normal CpG blocks\" visibility=full color=200,100,0 altColor=0,100,200 priority=20");
		tumorPvalueUnjustedWriter.println("track type=bedGraph name=\"-Log10 unjusted P value tumor CpG blocks \" description=\"-Log10 unjusted P value tumor CpG blocks\" visibility=full color=200,100,0 altColor=0,100,200 priority=20");
		String normalPvalueUnjustedCpGblockOverlapFn = prefix + "normalMerge_CpGblock_overlap.PvalueUnjusted." + minAlleleCpgReads + ".bedGraph";
		String tumorPvalueUnjustedCpGblockOverlapFn = prefix + "tumorMerge_CpGblock_overlap.PvalueUnjusted." + minAlleleCpgReads + ".bedGraph";
		PrintWriter normalPvalueUnjustedCpGblockOverlapWriter = new PrintWriter(new File(normalPvalueUnjustedCpGblockOverlapFn));
		PrintWriter tumorPvalueUnjustedCpGblockOverlapWriter = new PrintWriter(new File(tumorPvalueUnjustedCpGblockOverlapFn));
		normalPvalueUnjustedWriter.println("track type=bedGraph name=\"-Log10 unjusted P value normal CpG blocks (Overlap) \" description=\"-Log10 unjusted P value normal CpG blocks (Overlap)\" visibility=full color=200,100,0 altColor=0,100,200 priority=20");
		tumorPvalueUnjustedWriter.println("track type=bedGraph name=\"-Log10 unjusted P value tumor CpGblocks (Overlap) \" description=\"-Log10 unjusted P value tumor CpG blocks (Overlap)\" visibility=full color=200,100,0 altColor=0,100,200 priority=20");
	
		
		String normalMethyratioFn = prefix + "normalMerge_CpGblock.MethyRatio" + "." + minAlleleCpgReads + ".bedGraph";
		String tumorMethyratioFn = prefix + "tumorMerge_CpGblock.MethyRatio" + "." + minAlleleCpgReads + ".bedGraph";
		PrintWriter normalMethyratioWriter = new PrintWriter(new File(normalMethyratioFn));
		PrintWriter tumorMethyratioWriter = new PrintWriter(new File(tumorMethyratioFn));
		normalMethyratioWriter.println("track type=bedGraph name=\"absolute Log2 Methy ratio normal CpG block \" description=\"absolute Log2 Methy ratio normal CpG block\" visibility=full color=200,100,0 altColor=0,100,200 priority=20");
		tumorMethyratioWriter.println("track type=bedGraph name=\"absolute Log2 Methy ratio tumor CpG block \" description=\"absolute Log2 Methy ratio tumor CpG block\" visibility=full color=200,100,0 altColor=0,100,200 priority=20");
		String normalMethyratioCpGblockOverlapFn = prefix + "normalMerge_CpGblock_overlap.MethyRatio" + "." + minAlleleCpgReads + ".bedGraph";
		String tumorMethyratioCpGblockOverlapFn = prefix + "tumorMerge_CpGblock_overlap.MethyRatio" + "." + minAlleleCpgReads + ".bedGraph";
		PrintWriter normalMethyratioCpGblockOverlapWriter = new PrintWriter(new File(normalMethyratioCpGblockOverlapFn));
		PrintWriter tumorMethyratioCpGblockOverlapWriter = new PrintWriter(new File(tumorMethyratioCpGblockOverlapFn));
		normalMethyratioCpGblockOverlapWriter.println("track type=bedGraph name=\"absolute Log2 Methy ratio normal CpG block (Overlap) \" description=\"absolute Log2 Methy ratio normal CpG block\" visibility=full color=200,100,0 altColor=0,100,200 priority=20");
		tumorMethyratioCpGblockOverlapWriter.println("track type=bedGraph name=\"absolute Log2 Methy ratio tumor CpG block (Overlap) \" description=\"absolute Log2 Methy ratio tumor CpG block\" visibility=full color=200,100,0 altColor=0,100,200 priority=20");

		String normalMethyAlleleAFn = prefix + "normalMerge_CpGblock.MethyAlleleA" + "." + minAlleleCpgReads + ".bedGraph";
		String tumorMethyAlleleAFn = prefix + "tumorMerge_CpGblock.MethyAlleleA" + "." + minAlleleCpgReads + ".bedGraph";
		PrintWriter normalMethyAlleleAWriter = new PrintWriter(new File(normalMethyAlleleAFn));
		PrintWriter tumorMethyAlleleAWriter = new PrintWriter(new File(tumorMethyAlleleAFn));
		normalMethyAlleleAWriter.println("track type=bedGraph name=\"Allele A Methylation value normal CpG block \" description=\"Allele A Methylation value normal CpG block\" visibility=full color=200,100,0 altColor=0,100,200 priority=20");
		tumorMethyAlleleAWriter.println("track type=bedGraph name=\"Allele A Methylation value tumor CpG block \" description=\"Allele A Methylation value tumor CpG block\" visibility=full color=200,100,0 altColor=0,100,200 priority=20");
		String normalMethyAlleleACpGblockOverlapFn = prefix + "normalMerge_CpGblock_overlap.MethyAlleleA" + "." + minAlleleCpgReads + ".bedGraph";
		String tumorMethyAlleleACpGblockOverlapFn = prefix + "tumorMerge_CpGblock_overlap.MethyAlleleA" + "." + minAlleleCpgReads + ".bedGraph";
		PrintWriter normalMethyAlleleACpGblockOverlapWriter = new PrintWriter(new File(normalMethyAlleleACpGblockOverlapFn));
		PrintWriter tumorMethyAlleleACpGblockOverlapWriter = new PrintWriter(new File(tumorMethyAlleleACpGblockOverlapFn));
		normalMethyAlleleACpGblockOverlapWriter.println("track type=bedGraph name=\"Allele A Methylation value normal CpG block (Overlap) \" description=\"Allele A Methylation value normal CpG block\" visibility=full color=200,100,0 altColor=0,100,200 priority=20");
		tumorMethyAlleleACpGblockOverlapWriter.println("track type=bedGraph name=\"Allele A Methylation value tumor CpG block (Overlap) \" description=\"Allele A Methylation value tumor CpG block\" visibility=full color=200,100,0 altColor=0,100,200 priority=20");
		String normalMethyAlleleBFn = prefix + "normalMerge_CpGblock.MethyAlleleB" + "." + minAlleleCpgReads + ".bedGraph";
		String tumorMethyAlleleBFn = prefix + "tumorMerge_CpGblock.MethyAlleleB" + "." + minAlleleCpgReads + ".bedGraph";
		PrintWriter normalMethyAlleleBWriter = new PrintWriter(new File(normalMethyAlleleBFn));
		PrintWriter tumorMethyAlleleBWriter = new PrintWriter(new File(tumorMethyAlleleBFn));
		normalMethyAlleleBWriter.println("track type=bedGraph name=\"Allele B Methylation value normal CpG block \" description=\"Allele B Methylation value normal CpG block\" visibility=full color=200,100,0 altColor=0,100,200 priority=20");
		tumorMethyAlleleBWriter.println("track type=bedGraph name=\"Allele B Methylation value tumor CpG block \" description=\"Allele B Methylation value tumor CpG block\" visibility=full color=200,100,0 altColor=0,100,200 priority=20");
		String normalMethyAlleleBCpGblockOverlapFn = prefix + "normalMerge_CpGblock_overlap.MethyAlleleB" + "." + minAlleleCpgReads + ".bedGraph";
		String tumorMethyAlleleBCpGblockOverlapFn = prefix + "tumorMerge_CpGblock_overlap.MethyAlleleB" + "." + minAlleleCpgReads + ".bedGraph";
		PrintWriter normalMethyAlleleBCpGblockOverlapWriter = new PrintWriter(new File(normalMethyAlleleBCpGblockOverlapFn));
		PrintWriter tumorMethyAlleleBCpGblockOverlapWriter = new PrintWriter(new File(tumorMethyAlleleBCpGblockOverlapFn));
		normalMethyAlleleBCpGblockOverlapWriter.println("track type=bedGraph name=\"Allele B Methylation value normal CpG block (Overlap) \" description=\"Allele B Methylation value normal CpG block\" visibility=full color=200,100,0 altColor=0,100,200 priority=20");
		tumorMethyAlleleBCpGblockOverlapWriter.println("track type=bedGraph name=\"Allele B Methylation value tumor CpG block (Overlap) \" description=\"Allele B Methylation value tumor CpG block\" visibility=full color=200,100,0 altColor=0,100,200 priority=20");

		
		String normalSigGtfFn = prefix + "normalMerge_sigASM_CpGblock.PvalueAdjusted." + pValue + ".MethyRatio." + methyRatio + "." + minAlleleCpgReads + ".gtf";
		String tumorSigGtfFn = prefix + "tumorMerge_sigASM_CpGblock.PvalueAdjusted." + pValue + ".MethyRatio." + methyRatio + "." + minAlleleCpgReads + ".gtf";
		PrintWriter normalSigGtfWriter = new PrintWriter(new File(normalSigGtfFn));
		PrintWriter tumorSigGtfWriter = new PrintWriter(new File(tumorSigGtfFn));
		String normalSigCpGblockOverlapGtfFn = prefix + "normalMerge_sigASM_CpGblock_overlap.PvalueAdjusted." + pValue + ".MethyRatio." + methyRatio + "." + minAlleleCpgReads + ".gtf";
		String tumorSigCpGblockOverlapGtfFn = prefix + "tumorMerge_sigASM_CpGblock_overlap.PvalueAdjusted." + pValue + ".MethyRatio." + methyRatio + "." + minAlleleCpgReads + ".gtf";
		PrintWriter normalSigCpGblockOverlapGtfWriter = new PrintWriter(new File(normalSigCpGblockOverlapGtfFn));
		PrintWriter tumorSigCpGblockOverlapGtfWriter = new PrintWriter(new File(tumorSigCpGblockOverlapGtfFn));
		String normalOnlySigCpGblockOverlapGtfFn = prefix + "normalMerge_sigASM_normalOnly_CpGblock_overlap.PvalueAdjusted." + pValue + ".MethyRatio." + methyRatio + ".AlterPvalue." + alterpValue + ".AlterMethyRatio." + alterMethyRatio + "." + minAlleleCpgReads + ".gtf";
		String tumorOnlySigCpGblockOverlapGtfFn = prefix + "tumorMerge_sigASM_tumorOnly_CpGblock_overlap.PvalueAdjusted." + pValue + ".MethyRatio." + methyRatio + ".AlterPvalue." + alterpValue + ".AlterMethyRatio." + alterMethyRatio + "." + minAlleleCpgReads + ".gtf";
		PrintWriter normalOnlySigCpGblockOverlapGtfWriter = new PrintWriter(new File(normalOnlySigCpGblockOverlapGtfFn));
		PrintWriter tumorOnlySigCpGblockOverlapGtfWriter = new PrintWriter(new File(tumorOnlySigCpGblockOverlapGtfFn));
		String normalTumorBothSigCpGblockOverlapGtfFn = prefix + "normal-tumorMerge_sigASM_NormalTumorBoth_CpGblock_overlap.PvalueAdjusted." + pValue + ".MethyRatio." + methyRatio + ".AlterPvalue." + alterpValue + ".AlterMethyRatio." + alterMethyRatio + "." + minAlleleCpgReads + ".gtf";
		PrintWriter normalTumorBothSigCpGblockOverlapGtfWriter = new PrintWriter(new File(normalTumorBothSigCpGblockOverlapGtfFn));

		
		//int count = 0;
		int normalCount = 0;
		int tumorCount = 0;
		int cpGblockOverlapCount = 0;
		int sigNormalCount = 0;
		int sigTumorCount = 0;
		int sigNormalCpGblockOverlapCount = 0;
		int sigTumorCpGblockOverlapCount = 0;
		int sigNormalOnlyCpGblockOverlapCount = 0;
		int sigTumorOnlyCpGblockOverlapCount = 0;
		int sigNormalTumorBothCpGblockOverlapCount = 0;
		
		
		for (final String chr : chrs){
			String normalSnpFn = prefix + "normalMerge_" + chr + "_pValue_CpGblock.txt";
			String tumorSnpFn = prefix + "tumorMerge_" + chr + "_pValue_CpGblock.txt";
			File normalFile = new File(snpFilePath, normalSnpFn);
			File tumorFile = new File(snpFilePath, tumorSnpFn);
			BufferedReader normalBr = new BufferedReader(new FileReader(normalFile));
			BufferedReader tumorBr = new BufferedReader(new FileReader(tumorFile));
			TreeMap<Integer, String> normalSnpPosition = new TreeMap<Integer, String>();
			TreeMap<Integer, String> tumorSnpPosition = new TreeMap<Integer, String>();
			TreeMap<Integer, String> normalSigPosition = new TreeMap<Integer, String>();
			TreeMap<Integer, String> tumorSigPosition = new TreeMap<Integer, String>();
			String line;
			while( (line = normalBr.readLine()) != null){
				String[] tmpArray = line.split("\t");
				if(tmpArray[0].equalsIgnoreCase("NA"))
					continue;
				String[] range = tmpArray[1].split("~");
				int range1 = Integer.parseInt(range[0]);
				int range2 = Integer.parseInt(range[1]);
				
				if(Integer.parseInt(tmpArray[3]) >= minAlleleCpgReads && Integer.parseInt(tmpArray[4]) >= minAlleleCpgReads){
					normalSnpPosition.put(Integer.parseInt(tmpArray[0]), line);
					normalWriter.println(chr + "\t" + line);
					normalGtfWriter.printf("%s\tnormal\tCpG_Block\t%d\t%d\t0\t+\t.\n",chr,range1,range2);
					normalPvalueWriter.printf("%s\t%d\t%d\t%s\n",chr,range1,range2,tmpArray[13]);
					normalPvalueUnjustedWriter.printf("%s\t%d\t%d\t%s\n",chr,range1,range2,tmpArray[11]);
					if(tmpArray[9].equalsIgnoreCase("NA") || tmpArray[9].equalsIgnoreCase("-Inf") || tmpArray[9].equalsIgnoreCase("Inf") || tmpArray[9].equalsIgnoreCase("NAN"))
						normalMethyratioWriter.printf("%s\t%d\t%d\t%s\n",chr,range1,range2,"10");
					else
						normalMethyratioWriter.printf("%s\t%d\t%d\t%s\n",chr,range1,range2,tmpArray[9]);
					normalMethyAlleleAWriter.printf("%s\t%d\t%d\t%s\n",chr,range1,range2,tmpArray[7]);
					normalMethyAlleleBWriter.printf("%s\t%d\t%d\t%s\n",chr,range1,range2,tmpArray[8]);
					if((tmpArray[9].equalsIgnoreCase("NA") || tmpArray[9].equalsIgnoreCase("NAN")) && Double.parseDouble(tmpArray[12]) < pValue){
						normalSigGtfWriter.printf("%s\tnormal\tASM_CpG_Block\t%d\t%d\t0\t+\t.\n",chr,range1,range2);
						sigNormalCount++;
					}	
					else if((tmpArray[9].equalsIgnoreCase("NA") || tmpArray[9].equalsIgnoreCase("NAN")) && Double.parseDouble(tmpArray[12]) >= pValue){
						continue;
					}	
					else if((tmpArray[9].equalsIgnoreCase("-Inf") || tmpArray[9].equalsIgnoreCase("Inf")) && Double.parseDouble(tmpArray[12]) >= pValue){
						continue;
					}	
					else if((tmpArray[9].equalsIgnoreCase("-Inf") || tmpArray[9].equalsIgnoreCase("Inf")) && Double.parseDouble(tmpArray[12]) < pValue){
						normalSigGtfWriter.printf("%s\tnormal\tASM_CpG_Block\t%d\t%d\t0\t+\t.\n",chr,range1,range2);
						sigNormalCount++;
					}	
					else if(Double.parseDouble(tmpArray[9]) >= Math.log10(methyRatio)/Math.log10(2) && Double.parseDouble(tmpArray[12]) < pValue){
						normalSigGtfWriter.printf("%s\tnormal\tASM_CpG_Block\t%d\t%d\t0\t+\t.\n",chr,range1,range2);
						sigNormalCount++;
					}
					
				}
				
			}
			normalBr.close();
			while( (line = tumorBr.readLine()) != null){
				String[] tmpArray = line.split("\t");
				if(tmpArray[0].equalsIgnoreCase("NA"))
					continue;
				String[] range = tmpArray[1].split("~");
				int range1 = Integer.parseInt(range[0]);
				int range2 = Integer.parseInt(range[1]);
				
				if(Integer.parseInt(tmpArray[3]) >= minAlleleCpgReads && Integer.parseInt(tmpArray[4]) >= minAlleleCpgReads){
					tumorSnpPosition.put(Integer.parseInt(tmpArray[0]), line);
					tumorWriter.println(chr + "\t" + line);
					tumorGtfWriter.printf("%s\ttumor\tCpG_Block\t%d\t%d\t0\t+\t.\n",chr,range1,range2);
					tumorPvalueWriter.printf("%s\t%d\t%d\t%s\n",chr,range1,range2,tmpArray[13]);
					tumorPvalueUnjustedWriter.printf("%s\t%d\t%d\t%s\n",chr,range1,range2,tmpArray[11]);
					if(tmpArray[9].equalsIgnoreCase("NA") || tmpArray[9].equalsIgnoreCase("-Inf") || tmpArray[9].equalsIgnoreCase("Inf") || tmpArray[9].equalsIgnoreCase("NAN"))
						tumorMethyratioWriter.printf("%s\t%d\t%d\t%s\n",chr,range1,range2,"10");
					else
						tumorMethyratioWriter.printf("%s\t%d\t%d\t%s\n",chr,range1,range2,tmpArray[9]);
					tumorMethyAlleleAWriter.printf("%s\t%d\t%d\t%s\n",chr,range1,range2,tmpArray[7]);
					tumorMethyAlleleBWriter.printf("%s\t%d\t%d\t%s\n",chr,range1,range2,tmpArray[8]);
					if((tmpArray[9].equalsIgnoreCase("NA") || tmpArray[9].equalsIgnoreCase("NAN")) && Double.parseDouble(tmpArray[12]) < pValue){
						tumorSigGtfWriter.printf("%s\ttumor\tASM_CpG_Block\t%d\t%d\t0\t+\t.\n",chr,range1,range2);
						sigTumorCount++;
					}	
					else if((tmpArray[9].equalsIgnoreCase("NA") || tmpArray[9].equalsIgnoreCase("NAN")) && Double.parseDouble(tmpArray[12]) >= pValue){
						continue;
					}	
					else if((tmpArray[9].equalsIgnoreCase("-Inf") || tmpArray[9].equalsIgnoreCase("Inf")) && Double.parseDouble(tmpArray[12]) >= pValue){
						continue;
					}	
					else if((tmpArray[9].equalsIgnoreCase("-Inf") || tmpArray[9].equalsIgnoreCase("Inf")) && Double.parseDouble(tmpArray[12]) < pValue){
						tumorSigGtfWriter.printf("%s\ttumor\tASM_CpG_Block\t%d\t%d\t0\t+\t.\n",chr,range1,range2);
						sigTumorCount++;
					}	
					else if(Double.parseDouble(tmpArray[9]) >= Math.log10(methyRatio)/Math.log10(2) && Double.parseDouble(tmpArray[12]) < pValue){
						tumorSigGtfWriter.printf("%s\ttumor\tASM_CpG_Block\t%d\t%d\t0\t+\t.\n",chr,range1,range2);
						sigTumorCount++;
					}
					
				}
				
			}
			tumorBr.close();
			
			normalCount += normalSnpPosition.keySet().size();
			tumorCount += tumorSnpPosition.keySet().size();
			
			Iterator<Integer> normalIt = normalSnpPosition.keySet().iterator();
			while(normalIt.hasNext()){
				int tempSnp = normalIt.next();
				if(tumorSnpPosition.containsKey(tempSnp)){
					String value = normalSnpPosition.get(tempSnp);
					String[] tmpArray = value.split("\t");
					String[] range = tmpArray[1].split("~");
					int range1 = Integer.parseInt(range[0]);
					int range2 = Integer.parseInt(range[1]);
					if(tmpArray[0].equalsIgnoreCase("NA"))
						continue;
					cpGblockOverlapCount++;
					normalOverlapWriter.println(chr + "\t" + value);
					overlapGtfWriter.printf("%s\tnormal-tumor\tOverlap_CpG_Block\t%d\t%d\t0\t+\t.\n",chr,range1,range2);
					normalPvalueCpGblockOverlapWriter.printf("%s\t%d\t%d\t%s\n",chr,range1,range2,tmpArray[13]);
					normalPvalueUnjustedCpGblockOverlapWriter.printf("%s\t%d\t%d\t%s\n",chr,range1,range2,tmpArray[11]);
					if(tmpArray[9].equalsIgnoreCase("NA") || tmpArray[9].equalsIgnoreCase("-Inf") || tmpArray[9].equalsIgnoreCase("Inf") || tmpArray[9].equalsIgnoreCase("NAN"))
						normalMethyratioCpGblockOverlapWriter.printf("%s\t%d\t%d\t%s\n",chr,range1,range2,"10");
					else
						normalMethyratioCpGblockOverlapWriter.printf("%s\t%d\t%d\t%s\n",chr,range1,range2,tmpArray[9]);
					normalMethyAlleleACpGblockOverlapWriter.printf("%s\t%d\t%d\t%s\n",chr,range1,range2,tmpArray[7]);
					normalMethyAlleleBCpGblockOverlapWriter.printf("%s\t%d\t%d\t%s\n",chr,range1,range2,tmpArray[8]);
					//if (Double.parseDouble(tmpArray[12]) != 1)
						//System.err.println(tmpArray[12]);
					//System.err.println(tmpArray[9]);
					if((tmpArray[9].equalsIgnoreCase("NA") || tmpArray[9].equalsIgnoreCase("NAN")) && Double.parseDouble(tmpArray[12]) < pValue){
						normalSigCpGblockOverlapGtfWriter.printf("%s\tnormal\tASM_CpG_Block\t%d\t%d\t0\t+\t.\n",chr,range1,range2);
						normalSigPvalueCpGblockOverlapWriter.printf("%s\t%d\t%d\t%s\n",chr,range1,range2,tmpArray[13]);
						normalSigPosition.put(Integer.parseInt(tmpArray[0]), value);
					}	
					else if((tmpArray[9].equalsIgnoreCase("NA") || tmpArray[9].equalsIgnoreCase("NAN")) && Double.parseDouble(tmpArray[12]) >= pValue){
						continue;
					}	
					else if((tmpArray[9].equalsIgnoreCase("-Inf") || tmpArray[9].equalsIgnoreCase("Inf")) && Double.parseDouble(tmpArray[12]) >= pValue){
						continue;
					}	
					else if((tmpArray[9].equalsIgnoreCase("-Inf") || tmpArray[9].equalsIgnoreCase("Inf")) && Double.parseDouble(tmpArray[12]) < pValue){
						normalSigCpGblockOverlapGtfWriter.printf("%s\tnormal\tASM_CpG_Block\t%d\t%d\t0\t+\t.\n",chr,range1,range2);
						normalSigPvalueCpGblockOverlapWriter.printf("%s\t%d\t%d\t%s\n",chr,range1,range2,tmpArray[13]);
						normalSigPosition.put(Integer.parseInt(tmpArray[0]), value);
					}	
					else if(Double.parseDouble(tmpArray[9]) >= Math.log10(methyRatio)/Math.log10(2) && Double.parseDouble(tmpArray[12]) < pValue){
						normalSigCpGblockOverlapGtfWriter.printf("%s\tnormal\tASM_CpG_Block\t%d\t%d\t0\t+\t.\n",chr,range1,range2);
						normalSigPvalueCpGblockOverlapWriter.printf("%s\t%d\t%d\t%s\n",chr,range1,range2,tmpArray[13]);
						normalSigPosition.put(Integer.parseInt(tmpArray[0]), value);
					}					
				}
			}

			Iterator<Integer> tumorIt = tumorSnpPosition.keySet().iterator();
			while(tumorIt.hasNext()){
				int tempSnp = tumorIt.next();
				if(normalSnpPosition.containsKey(tempSnp)){
					String value = tumorSnpPosition.get(tempSnp);
					String[] tmpArray = value.split("\t");
					String[] range = tmpArray[1].split("~");
					int range1 = Integer.parseInt(range[0]);
					int range2 = Integer.parseInt(range[1]);
					if(tmpArray[0].equalsIgnoreCase("NA"))
						continue;
					tumorOverlapWriter.println(chr + "\t" + value);
					tumorPvalueCpGblockOverlapWriter.printf("%s\t%d\t%d\t%s\n",chr,range1,range2,tmpArray[13]);
					tumorPvalueUnjustedCpGblockOverlapWriter.printf("%s\t%d\t%d\t%s\n",chr,range1,range2,tmpArray[11]);
					if(tmpArray[9].equalsIgnoreCase("NA") || tmpArray[9].equalsIgnoreCase("-Inf") || tmpArray[9].equalsIgnoreCase("Inf") || tmpArray[9].equalsIgnoreCase("NAN"))
						tumorMethyratioCpGblockOverlapWriter.printf("%s\t%d\t%d\t%s\n",chr,range1,range2,"10");
					else
						tumorMethyratioCpGblockOverlapWriter.printf("%s\t%d\t%d\t%s\n",chr,range1,range2,tmpArray[9]);
					tumorMethyAlleleACpGblockOverlapWriter.printf("%s\t%d\t%d\t%s\n",chr,range1,range2,tmpArray[7]);
					tumorMethyAlleleBCpGblockOverlapWriter.printf("%s\t%d\t%d\t%s\n",chr,range1,range2,tmpArray[8]);
					if((tmpArray[9].equalsIgnoreCase("NA") || tmpArray[9].equalsIgnoreCase("NAN")) && Double.parseDouble(tmpArray[12]) < pValue){
						tumorSigCpGblockOverlapGtfWriter.printf("%s\ttumor\tASM_CpG_Block\t%d\t%d\t0\t+\t.\n",chr,range1,range2);
						tumorSigPvalueCpGblockOverlapWriter.printf("%s\t%d\t%d\t%s\n",chr,range1,range2,tmpArray[13]);
						tumorSigPosition.put(Integer.parseInt(tmpArray[0]), value);
					}	
					else if((tmpArray[9].equalsIgnoreCase("NA") || tmpArray[9].equalsIgnoreCase("NAN")) && Double.parseDouble(tmpArray[12]) >= pValue){
						continue;
					}	
					else if((tmpArray[9].equalsIgnoreCase("-Inf") || tmpArray[9].equalsIgnoreCase("Inf")) && Double.parseDouble(tmpArray[12]) >= pValue){
						continue;
					}	
					else if((tmpArray[9].equalsIgnoreCase("-Inf") || tmpArray[9].equalsIgnoreCase("Inf")) && Double.parseDouble(tmpArray[12]) < pValue){
						tumorSigCpGblockOverlapGtfWriter.printf("%s\ttumor\tASM_CpG_Block\t%d\t%d\t0\t+\t.\n",chr,range1,range2);
						tumorSigPvalueCpGblockOverlapWriter.printf("%s\t%d\t%d\t%s\n",chr,range1,range2,tmpArray[13]);
						tumorSigPosition.put(Integer.parseInt(tmpArray[0]), value);
					}	
					else if(Double.parseDouble(tmpArray[9]) >= Math.log10(methyRatio)/Math.log10(2) && Double.parseDouble(tmpArray[12]) < pValue){
						tumorSigCpGblockOverlapGtfWriter.printf("%s\ttumor\tASM_CpG_Block\t%d\t%d\t0\t+\t.\n",chr,range1,range2);
						tumorSigPvalueCpGblockOverlapWriter.printf("%s\t%d\t%d\t%s\n",chr,range1,range2,tmpArray[13]);
						tumorSigPosition.put(Integer.parseInt(tmpArray[0]), value);
					}
				}
			}
			
			sigNormalCpGblockOverlapCount += normalSigPosition.keySet().size();
			sigTumorCpGblockOverlapCount += tumorSigPosition.keySet().size();
			
			Iterator<Integer> normalSigIt = normalSigPosition.keySet().iterator();
			while(normalSigIt.hasNext()){
				int tempSnp = normalSigIt.next();
				String value = normalSnpPosition.get(tempSnp);
				String[] tmpArray = value.split("\t");
				String[] range = tmpArray[1].split("~");
				int range1 = Integer.parseInt(range[0]);
				int range2 = Integer.parseInt(range[1]);
				if(tumorSigPosition.containsKey(tempSnp)){
					normalTumorBothSigCpGblockOverlapGtfWriter.printf("%s\tnormal-tumor\tboth_ASM_CpG_Block\t%d\t%d\t0\t+\t.\n",chr,range1,range2);
					sigNormalTumorBothCpGblockOverlapCount++;
				}
				else{
					
					String tmpLine = tumorSnpPosition.get(tempSnp);
					String[] tmpControlArray = tmpLine.split("\t");
					if(((tmpArray[9].equalsIgnoreCase("-Inf") || tmpArray[9].equalsIgnoreCase("Inf") || tmpArray[9].equalsIgnoreCase("NA") || tmpArray[9].equalsIgnoreCase("NAN")) && (!tmpControlArray[9].equalsIgnoreCase("-Inf") && !tmpControlArray[9].equalsIgnoreCase("Inf") && !tmpControlArray[9].equalsIgnoreCase("NA") && !tmpControlArray[9].equalsIgnoreCase("NAN"))) && Double.parseDouble(tmpControlArray[10]) >= alterpValue ){
						normalOnlySigCpGblockOverlapGtfWriter.printf("%s\tnormal\tnormalOnly_ASM_CpG_Block\t%d\t%d\t0\t+\t.\n",chr,range1,range2);
						sigNormalOnlyCpGblockOverlapCount++;
					}
					else if(((tmpControlArray[9].equalsIgnoreCase("-Inf") || tmpControlArray[9].equalsIgnoreCase("Inf") || tmpControlArray[9].equalsIgnoreCase("NA")) && (!tmpArray[9].equalsIgnoreCase("-Inf") && !tmpArray[9].equalsIgnoreCase("Inf") && !tmpArray[9].equalsIgnoreCase("NA"))) && Double.parseDouble(tmpControlArray[10]) >= alterpValue ){
						normalOnlySigCpGblockOverlapGtfWriter.printf("%s\tnormal\tnormalOnly_ASM_CpG_Block\t%d\t%d\t0\t+\t.\n",chr,range1,range2);
						sigNormalOnlyCpGblockOverlapCount++;
					}
					else if(((tmpControlArray[9].equalsIgnoreCase("-Inf") || tmpControlArray[9].equalsIgnoreCase("Inf") || tmpControlArray[9].equalsIgnoreCase("NA")) || (tmpArray[9].equalsIgnoreCase("-Inf") || tmpArray[9].equalsIgnoreCase("Inf") || tmpArray[9].equalsIgnoreCase("NA")))){
						continue;
					}
					else if(Double.parseDouble(tmpArray[9]) <= Math.abs(Math.log10(1.5)/Math.log10(2)) && Double.parseDouble(tmpControlArray[6]) >= alterpValue){
						normalOnlySigCpGblockOverlapGtfWriter.printf("%s\tnormal\tnormalOnly_ASM_CpG_Block\t%d\t%d\t0\t+\t.\n",chr,range1,range2);
						sigNormalOnlyCpGblockOverlapCount++;
					}
				}
			}
			
			Iterator<Integer> tumorSigIt = tumorSigPosition.keySet().iterator();
			while(tumorSigIt.hasNext()){
				int tempSnp = tumorSigIt.next();
				String value = tumorSnpPosition.get(tempSnp);
				String[] tmpArray = value.split("\t");
				String[] range = tmpArray[1].split("~");
				int range1 = Integer.parseInt(range[0]);
				int range2 = Integer.parseInt(range[1]);
				if(normalSigPosition.containsKey(tempSnp)){
					
				}
				else{
					
					String tmpLine = normalSnpPosition.get(tempSnp);
					String[] tmpControlArray = tmpLine.split("\t");
					if(((tmpArray[9].equalsIgnoreCase("-Inf") || tmpArray[9].equalsIgnoreCase("Inf") || tmpArray[9].equalsIgnoreCase("NA") || tmpArray[9].equalsIgnoreCase("NAN")) && (!tmpControlArray[9].equalsIgnoreCase("-Inf") && !tmpControlArray[9].equalsIgnoreCase("Inf") && !tmpControlArray[9].equalsIgnoreCase("NA") && !tmpControlArray[9].equalsIgnoreCase("NAN"))) && Double.parseDouble(tmpControlArray[10]) >= alterpValue ){
						tumorOnlySigCpGblockOverlapGtfWriter.printf("%s\ttumor\ttumorOnly_ASM_CpG_Block\t%d\t%d\t0\t+\t.\n",chr,range1,range2);
						sigTumorOnlyCpGblockOverlapCount++;
					}
					else if(((tmpControlArray[9].equalsIgnoreCase("-Inf") || tmpControlArray[9].equalsIgnoreCase("Inf") || tmpControlArray[9].equalsIgnoreCase("NA")) && (!tmpArray[9].equalsIgnoreCase("-Inf") && !tmpArray[9].equalsIgnoreCase("Inf") && !tmpArray[9].equalsIgnoreCase("NA"))) && Double.parseDouble(tmpControlArray[10]) >= alterpValue ){
						tumorOnlySigCpGblockOverlapGtfWriter.printf("%s\ttumor\ttumorOnly_ASM_CpG_Block\t%d\t%d\t0\t+\t.\n",chr,range1,range2);
						sigTumorOnlyCpGblockOverlapCount++;
					}
					else if(((tmpControlArray[9].equalsIgnoreCase("-Inf") || tmpControlArray[9].equalsIgnoreCase("Inf") || tmpControlArray[9].equalsIgnoreCase("NA")) || (tmpArray[9].equalsIgnoreCase("-Inf") || tmpArray[9].equalsIgnoreCase("Inf") || tmpArray[9].equalsIgnoreCase("NA")))){
						continue;
					}
					else if(Double.parseDouble(tmpArray[9]) <= Math.abs(Math.log10(1.5)/Math.log10(2)) && Double.parseDouble(tmpControlArray[6]) >= alterpValue){
						tumorOnlySigCpGblockOverlapGtfWriter.printf("%s\ttumor\ttumorOnly_ASM_CpG_Block\t%d\t%d\t0\t+\t.\n",chr,range1,range2);
						sigTumorOnlyCpGblockOverlapCount++;
					}
				}
			}
		}

		normalWriter.close();
		tumorWriter.close();
		normalGtfWriter.close();
		tumorGtfWriter.close();
		normalOverlapWriter.close();
		tumorOverlapWriter.close();
		overlapGtfWriter.close();
		normalPvalueWriter.close();
		tumorPvalueWriter.close();
		normalPvalueCpGblockOverlapWriter.close();
		tumorPvalueCpGblockOverlapWriter.close();
		normalSigPvalueCpGblockOverlapWriter.close();
		tumorSigPvalueCpGblockOverlapWriter.close();
		normalPvalueUnjustedWriter.close();
		tumorPvalueUnjustedWriter.close();
		normalPvalueUnjustedCpGblockOverlapWriter.close();
		tumorPvalueUnjustedCpGblockOverlapWriter.close();
		normalMethyratioWriter.close();
		tumorMethyratioWriter.close();
		normalMethyratioCpGblockOverlapWriter.close();
		tumorMethyratioCpGblockOverlapWriter.close();
		normalMethyAlleleAWriter.close();
		tumorMethyAlleleAWriter.close();
		normalMethyAlleleACpGblockOverlapWriter.close();
		tumorMethyAlleleACpGblockOverlapWriter.close();
		normalMethyAlleleBWriter.close();
		tumorMethyAlleleBWriter.close();
		normalMethyAlleleBCpGblockOverlapWriter.close();
		tumorMethyAlleleBCpGblockOverlapWriter.close();
		normalSigGtfWriter.close();
		tumorSigGtfWriter.close();
		normalSigCpGblockOverlapGtfWriter.close();
		tumorSigCpGblockOverlapGtfWriter.close();
		normalOnlySigCpGblockOverlapGtfWriter.close();
		tumorOnlySigCpGblockOverlapGtfWriter.close();
		normalTumorBothSigCpGblockOverlapGtfWriter.close();
		
		
		System.err.printf("The count of normal CpG Block is: %d\n",normalCount);
		System.err.printf("The count of tumor CpG Block is: %d\n",tumorCount);
		System.err.printf("The count of overlap CpG block is: %d\n",cpGblockOverlapCount);
		System.err.printf("The count of significant normal ASM Block is: %d\n", sigNormalCount);
		System.err.printf("The count of significant tumor ASM Block is: %d\n", sigTumorCount);
		System.err.printf("The count of significant normal ASM Block in overlap CpG block is: %d\n", sigNormalCpGblockOverlapCount);
		System.err.printf("The count of significant tumor ASM Block in overlap CpG block is: %d\n", sigTumorCpGblockOverlapCount);
		System.err.printf("The count of both significant ASM block in overlap CpG block is: %d\n",sigNormalTumorBothCpGblockOverlapCount);
		System.err.printf("The count of significant changed normal ASM Block in overlap CpG block is: %d\n", sigNormalOnlyCpGblockOverlapCount);
		System.err.printf("The count of significant changed tumor ASM Block in overlap CpG block is: %d\n", sigTumorOnlyCpGblockOverlapCount);
	}

}
