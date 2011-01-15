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

public class OverlapAsmBlock {

	/**
	 * @param args
	 */

	final private static String USAGE = "OverlapSnp [opt] ASMBlockFileDirecotryPath";
	final private static String prefix = "methylCGsRich_ASM_AllSnp_";
	
	@Option(name="-chrom",multiValued=true,usage="One or more chroms, eg. --chrom chr1 --chrom chr5")
	protected List<String> chrs = new ArrayList<String>(25);
	@Option(name="-cpgNumInBlock",usage="minimum number in one ASM Block, default value is: 5")
	protected int cpgNumInBlock = 5;
	@Option(name="-pValue",usage="p value thresh hold for significant ASM block, default value is: 0.01")
	protected double pValue = 0.01;
	@Option(name="-methyRatio",usage="methylation ratio thresh hold for significant ASM block, default value is: 2")
	protected double methyRatio = 2;
	@Option(name="-alterpValue",usage="p value thresh hold for not significant ASM block, default value is: 0.05")
	protected double alterpValue = 0.05;
	@Option(name="-alterMethyRatio",usage="methylation ratio thresh hold for not significant ASM block, default value is: 1.5")
	protected double alterMethyRatio = 1.5;
	
	@Argument
	private List<String> stringArgs = new ArrayList<String>();
	
	public static void main(String[] args)
	throws Exception
	{
		new OverlapAsmBlock().doMain(args);
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


		String normalFn = prefix + "normalMerge_pValue_ASMblock_overlap." + cpgNumInBlock + ".txt";
		String tumorFn = prefix + "tumorMerge_pValue_ASMblock_overlap." + cpgNumInBlock + ".txt";
		PrintWriter normalWriter = new PrintWriter(new File(normalFn));
		PrintWriter tumorWriter = new PrintWriter(new File(tumorFn));
		
		String normalPvalueFn = prefix + "normalMerge_ASMblock_overlap.PvalueAdjusted." + cpgNumInBlock + ".bedGraph";
		String tumorPvalueFn = prefix + "tumorMerge_ASMblock_overlap.PvalueAdjusted." + cpgNumInBlock + ".bedGraph";
		PrintWriter normalPvalueWriter = new PrintWriter(new File(normalPvalueFn));
		PrintWriter tumorPvalueWriter = new PrintWriter(new File(tumorPvalueFn));
		normalPvalueWriter.println("track type=bedGraph name=\"-Log10 P value normal ASM block \" description=\"-Log10 P value normal ASM block\" visibility=full color=200,100,0 altColor=0,100,200 priority=20");
		tumorPvalueWriter.println("track type=bedGraph name=\"-Log10 P value tumor ASM block \" description=\"-Log10 P value tumor ASM block\" visibility=full color=200,100,0 altColor=0,100,200 priority=20");
		
		
		String normalPvalueUnjustFn = prefix + "normalMerge_ASMblock_overlap.PvalueUnjust" + "." + cpgNumInBlock + ".bedGraph";
		String tumorPvalueUnjustFn = prefix + "tumorMerge_ASMblock_overlap.PvalueUnjust" +  "." + cpgNumInBlock + ".bedGraph";
		PrintWriter normalPvalueUnjustWriter = new PrintWriter(new File(normalPvalueUnjustFn));
		PrintWriter tumorPvalueUnjustWriter = new PrintWriter(new File(tumorPvalueUnjustFn));
		normalPvalueUnjustWriter.println("track type=bedGraph name=\"-Log10 P value Unjusted normal ASM block \" description=\"-Log10 P value Unjusted normal ASM block\" visibility=full color=200,100,0 altColor=0,100,200 priority=20");
		tumorPvalueUnjustWriter.println("track type=bedGraph name=\"-Log10 P value Unjusted tumor ASM block \" description=\"-Log10 P value Unjusted tumor ASM block\" visibility=full color=200,100,0 altColor=0,100,200 priority=20");

		
		String normalMethyFn = prefix + "normalMerge_ASMblock_overlap.Methy" + "." + cpgNumInBlock + ".bedGraph";
		String tumorMethyFn = prefix + "tumorMerge_ASMblock_overlap.Methy" + "." + cpgNumInBlock + ".bedGraph";
		PrintWriter normalMethyWriter = new PrintWriter(new File(normalMethyFn));
		PrintWriter tumorMethyWriter = new PrintWriter(new File(tumorMethyFn));
		normalMethyWriter.println("track type=bedGraph name=\"Log2 Methy ratio normal ASM block \" description=\"Log2 Methy ratio normal ASM block\" visibility=full color=200,100,0 altColor=0,100,200 priority=20");
		tumorMethyWriter.println("track type=bedGraph name=\"Log2 Methy ratio tumor ASM block \" description=\"Log2 Methy ratio tumor ASM block\" visibility=full color=200,100,0 altColor=0,100,200 priority=20");

		String normalSigGtfFn = prefix + "normalMerge_ASMblock_overlap.sigGroup.AdjustedPvalue." + pValue + "." + methyRatio + "." + cpgNumInBlock + ".gtf";
		String tumorSigGtfFn = prefix + "tumorMerge_ASMblock_overlap.sigGroup.AdjustedPvalue." + pValue + "." + methyRatio + "." + cpgNumInBlock + ".gtf";
		PrintWriter normalSigGtfWriter = new PrintWriter(new File(normalSigGtfFn));
		PrintWriter tumorSigGtfWriter = new PrintWriter(new File(tumorSigGtfFn));
		
		String normalSigGtfNonoverlapFn = prefix + "normalMerge_ASMblock_overlap.sigGroup.Nonoverlap.AdjustedPvalue." + pValue + "." + methyRatio + "." + cpgNumInBlock + ".gtf";
		String tumorSigGtfNonoverlapFn = prefix + "tumorMerge_ASMblock_overlap.sigGroup.Nonoverlap.AdjustedPvalue." + pValue + "." + methyRatio + "." + cpgNumInBlock + ".gtf";
		String sigGtfOverlapFn = prefix + "normal-tumorMerge_ASMblock_overlap.sigGroup.Overlap.AdjustedPvalue." + pValue + "." + methyRatio + "." + cpgNumInBlock + ".gtf";
		PrintWriter normalSigGtfNonoverlapWriter = new PrintWriter(new File(normalSigGtfNonoverlapFn));
		PrintWriter tumorSigGtfNonoverlapWriter = new PrintWriter(new File(tumorSigGtfNonoverlapFn));
		PrintWriter sigGtfOverlapWriter = new PrintWriter(new File(sigGtfOverlapFn));
		
		String normalSigGtfChangeFn = prefix + "normalMerge_ASMblock_overlap.sigGroup.sigChange(1,0.01).AdjustedPvalue." + pValue + "." + methyRatio + "." + cpgNumInBlock + ".gtf";
		String tumorSigGtfChangeFn = prefix + "tumorMerge_ASMblock_overlap.sigGroup.sigChange(1,0.01).AdjustedPvalue." + pValue + "." + methyRatio + "." + cpgNumInBlock + ".gtf";
		PrintWriter normalSigGtfChangeWriter = new PrintWriter(new File(normalSigGtfChangeFn));
		PrintWriter tumorSigGtfChangeWriter = new PrintWriter(new File(tumorSigGtfChangeFn));
		
		int count = 0;
		int normalCount = 0;
		int tumorCount = 0;
		int sigNormalCount = 0;
		int sigTumorCount = 0;
		int sigOverlapCount = 0;
		int sigNormalChangeCount = 0;
		int sigTumorChangeCount = 0;
		
		for (final String chr : chrs){
			String normalSnpFn = prefix + "normalMerge_" + chr + "_pValue_ASMblock.txt";
			String tumorSnpFn = prefix + "tumorMerge_" + chr + "_pValue_ASMblock.txt";
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
				if(Integer.parseInt(tmpArray[2]) >= cpgNumInBlock){
					normalSnpPosition.put(Integer.parseInt(tmpArray[0]), line);
				}
				
			}
			normalBr.close();
			while( (line = tumorBr.readLine()) != null){
				String[] tmpArray = line.split("\t");
				if(tmpArray[0].equalsIgnoreCase("NA"))
					continue;
				if(Integer.parseInt(tmpArray[2]) >= cpgNumInBlock){
					tumorSnpPosition.put(Integer.parseInt(tmpArray[0]), line);
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
					normalWriter.println(chr + "\t" + value);
					count++;
					normalPvalueWriter.printf("%s\t%d\t%d\t%s\n",chr,range1,range2,tmpArray[7]);
					if(tmpArray[3].equalsIgnoreCase("NA") || tmpArray[3].equalsIgnoreCase("-Inf") || tmpArray[3].equalsIgnoreCase("Inf"))
						normalMethyWriter.printf("%s\t%d\t%d\t%s\n",chr,range1,range2,"10");
					else
						normalMethyWriter.printf("%s\t%d\t%d\t%s\n",chr,range1,range2,tmpArray[3]);
					normalPvalueUnjustWriter.printf("%s\t%d\t%d\t%s\n",chr,range1,range2,tmpArray[5]);
					if(tmpArray[3].equalsIgnoreCase("NA") && Double.parseDouble(tmpArray[6]) < pValue){
						normalSigGtfWriter.printf("%s\tnormal\tASM_Block\t%d\t%d\t0\t+\t.\n",chr,range1,range2);
						normalSigPosition.put(Integer.parseInt(tmpArray[0]), value);
					}	
					else if(tmpArray[3].equalsIgnoreCase("NA") && Double.parseDouble(tmpArray[6]) >= pValue){
						continue;
					}	
					else if((tmpArray[3].equalsIgnoreCase("-Inf") || tmpArray[3].equalsIgnoreCase("Inf")) && Double.parseDouble(tmpArray[6]) >= pValue){
						continue;
					}	
					else if((tmpArray[3].equalsIgnoreCase("-Inf") || tmpArray[3].equalsIgnoreCase("Inf")) && Double.parseDouble(tmpArray[6]) < pValue){
						normalSigGtfWriter.printf("%s\tnormal\tASM_Block\t%d\t%d\t0\t+\t.\n",chr,range1,range2);
						normalSigPosition.put(Integer.parseInt(tmpArray[0]), value);
					}	
					else if((Double.parseDouble(tmpArray[3]) >= Math.log10(methyRatio)/Math.log10(2) || Double.parseDouble(tmpArray[3]) <= -Math.log10(methyRatio)/Math.log10(2)) && Double.parseDouble(tmpArray[6]) < pValue){
						normalSigGtfWriter.printf("%s\tnormal\tASM_Block\t%d\t%d\t0\t+\t.\n",chr,range1,range2);
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
					tumorWriter.println(chr + "\t" + value);
					tumorPvalueWriter.printf("%s\t%d\t%d\t%s\n",chr,range1,range2,tmpArray[7]);
					if(tmpArray[3].equalsIgnoreCase("NA") || tmpArray[3].equalsIgnoreCase("-Inf") || tmpArray[3].equalsIgnoreCase("Inf"))
						tumorMethyWriter.printf("%s\t%d\t%d\t%s\n",chr,range1,range2,"10");
					else
						tumorMethyWriter.printf("%s\t%d\t%d\t%s\n",chr,range1,range2,tmpArray[3]);
					tumorPvalueUnjustWriter.printf("%s\t%d\t%d\t%s\n",chr,range1,range2,tmpArray[5]);
					if(tmpArray[3].equalsIgnoreCase("NA") && Double.parseDouble(tmpArray[6]) < pValue){
						tumorSigGtfWriter.printf("%s\ttumor\tASM_Block\t%d\t%d\t0\t+\t.\n",chr,range1,range2);
						tumorSigPosition.put(Integer.parseInt(tmpArray[0]), value);
					}	
					else if(tmpArray[3].equalsIgnoreCase("NA") && Double.parseDouble(tmpArray[6]) >= pValue){
						continue;
					}
					else if((tmpArray[3].equalsIgnoreCase("-Inf") || tmpArray[3].equalsIgnoreCase("Inf")) && Double.parseDouble(tmpArray[6]) >= pValue){
						continue;
					}	
					else if((tmpArray[3].equalsIgnoreCase("-Inf") || tmpArray[3].equalsIgnoreCase("Inf")) && Double.parseDouble(tmpArray[6]) < pValue){
						tumorSigGtfWriter.printf("%s\ttumor\tASM_Block\t%d\t%d\t0\t+\t.\n",chr,range1,range2);
						tumorSigPosition.put(Integer.parseInt(tmpArray[0]), value);
					}	
					else if((Double.parseDouble(tmpArray[3]) >= Math.log10(methyRatio)/Math.log10(methyRatio) || Double.parseDouble(tmpArray[3]) <= -Math.log10(methyRatio)/Math.log10(methyRatio)) && Double.parseDouble(tmpArray[6]) < pValue){
						tumorSigGtfWriter.printf("%s\ttumor\tASM_Block\t%d\t%d\t0\t+\t.\n",chr,range1,range2);
						tumorSigPosition.put(Integer.parseInt(tmpArray[0]), value);
					}
				}
			}
			
			sigNormalCount += normalSigPosition.keySet().size();
			sigTumorCount += tumorSigPosition.keySet().size();
			
			Iterator<Integer> normalSigIt = normalSigPosition.keySet().iterator();
			while(normalSigIt.hasNext()){
				int tempSnp = normalSigIt.next();
				String value = normalSnpPosition.get(tempSnp);
				String[] tmpArray = value.split("\t");
				String[] range = tmpArray[1].split("~");
				int range1 = Integer.parseInt(range[0]);
				int range2 = Integer.parseInt(range[1]);
				if(tumorSigPosition.containsKey(tempSnp)){
					sigGtfOverlapWriter.printf("%s\tnormal-tumor\tASM_Block_overlap\t%d\t%d\t0\t+\t.\n",chr,range1,range2);
					sigOverlapCount++;
				}
				else{
					normalSigGtfNonoverlapWriter.printf("%s\tnormal\tASM_Block_nonoverlap\t%d\t%d\t0\t+\t.\n",chr,range1,range2);
					String tmpLine = tumorSnpPosition.get(tempSnp);
					String[] tmpControlArray = tmpLine.split("\t");
					if(((tmpArray[3].equalsIgnoreCase("-Inf") || tmpArray[3].equalsIgnoreCase("Inf") || tmpArray[3].equalsIgnoreCase("NA")) && (!tmpControlArray[3].equalsIgnoreCase("-Inf") && !tmpControlArray[3].equalsIgnoreCase("Inf") && !tmpControlArray[3].equalsIgnoreCase("NA"))) && Double.parseDouble(tmpControlArray[6]) >= 0.01 ){
						normalSigGtfChangeWriter.printf("%s\tnormal\tASM_Block_dynamicChange\t%d\t%d\t0\t+\t.\n",chr,range1,range2);
						sigNormalChangeCount++;
					}
					else if(((tmpControlArray[3].equalsIgnoreCase("-Inf") || tmpControlArray[3].equalsIgnoreCase("Inf") || tmpControlArray[3].equalsIgnoreCase("NA")) && (!tmpArray[3].equalsIgnoreCase("-Inf") && !tmpArray[3].equalsIgnoreCase("Inf") && !tmpArray[3].equalsIgnoreCase("NA"))) && Double.parseDouble(tmpControlArray[6]) >= 0.01 ){
						normalSigGtfChangeWriter.printf("%s\tnormal\tASM_Block_dynamicChange\t%d\t%d\t0\t+\t.\n",chr,range1,range2);
						sigNormalChangeCount++;
					}
					else if(((tmpControlArray[3].equalsIgnoreCase("-Inf") || tmpControlArray[3].equalsIgnoreCase("Inf") || tmpControlArray[3].equalsIgnoreCase("NA")) || (tmpArray[3].equalsIgnoreCase("-Inf") || tmpArray[3].equalsIgnoreCase("Inf") || tmpArray[3].equalsIgnoreCase("NA")))){
						continue;
					}
					else if(Math.abs(Math.abs(Double.parseDouble(tmpArray[3])) - Math.abs(Double.parseDouble(tmpControlArray[3]))) >= 1 && Double.parseDouble(tmpControlArray[6]) >= 0.01){
						normalSigGtfChangeWriter.printf("%s\tnormal\tASM_Block_dynamicChange\t%d\t%d\t0\t+\t.\n",chr,range1,range2);
						sigNormalChangeCount++;
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
					tumorSigGtfNonoverlapWriter.printf("%s\ttumor\tASM_Block_nonoverlap\t%d\t%d\t0\t+\t.\n",chr,range1,range2);
					String tmpLine = normalSnpPosition.get(tempSnp);
					String[] tmpControlArray = tmpLine.split("\t");
					if(((tmpArray[3].equalsIgnoreCase("-Inf") || tmpArray[3].equalsIgnoreCase("Inf") || tmpArray[3].equalsIgnoreCase("NA")) && (!tmpControlArray[3].equalsIgnoreCase("-Inf") && !tmpControlArray[3].equalsIgnoreCase("Inf") && !tmpControlArray[3].equalsIgnoreCase("NA"))) && Double.parseDouble(tmpControlArray[6]) >= 0.01 ){
						tumorSigGtfChangeWriter.printf("%s\ttumor\tASM_Block_dynamicChange\t%d\t%d\t0\t+\t.\n",chr,range1,range2);
						sigTumorChangeCount++;
					}
					else if(((tmpControlArray[3].equalsIgnoreCase("-Inf") || tmpControlArray[3].equalsIgnoreCase("Inf") || tmpControlArray[3].equalsIgnoreCase("NA")) && (!tmpArray[3].equalsIgnoreCase("-Inf") && !tmpArray[3].equalsIgnoreCase("Inf") && !tmpArray[3].equalsIgnoreCase("NA"))) && Double.parseDouble(tmpControlArray[6]) >= 0.01 ){
						tumorSigGtfChangeWriter.printf("%s\ttumor\tASM_Block_dynamicChange\t%d\t%d\t0\t+\t.\n",chr,range1,range2);
						sigTumorChangeCount++;
					}
					else if(((tmpControlArray[3].equalsIgnoreCase("-Inf") || tmpControlArray[3].equalsIgnoreCase("Inf") || tmpControlArray[3].equalsIgnoreCase("NA")) || (tmpArray[3].equalsIgnoreCase("-Inf") || tmpArray[3].equalsIgnoreCase("Inf") || tmpArray[3].equalsIgnoreCase("NA")))){
						continue;
					}
					else if(Math.abs(Math.abs(Double.parseDouble(tmpArray[3])) - Math.abs(Double.parseDouble(tmpControlArray[3]))) >= 1 && Double.parseDouble(tmpControlArray[6]) >= 0.01){
						tumorSigGtfChangeWriter.printf("%s\ttumor\tASM_Block_dynamicChange\t%d\t%d\t0\t+\t.\n",chr,range1,range2);
						sigTumorChangeCount++;
					}
				}
			}
		}

		normalWriter.close();
		normalPvalueWriter.close();
		normalMethyWriter.close();
		normalSigGtfWriter.close();
		normalPvalueUnjustWriter.close();
		//Iterator<Integer> it2 = snpPosition2.iterator();
		
		tumorWriter.close();
		tumorPvalueWriter.close();
		tumorMethyWriter.close();
		tumorSigGtfWriter.close();
		tumorPvalueUnjustWriter.close();
		
		sigGtfOverlapWriter.close();
		normalSigGtfNonoverlapWriter.close();
		tumorSigGtfNonoverlapWriter.close();
		normalSigGtfChangeWriter.close();
		tumorSigGtfChangeWriter.close();
		
		System.err.printf("The count of normal ASM Block (contains %d CpGs) is: %d\n",cpgNumInBlock, normalCount);
		System.err.printf("The count of tumor ASM Block (contains %d CpGs) is: %d\n",cpgNumInBlock, tumorCount);
		System.err.printf("The count of overlap ASM block is: %d\n",count);
		System.err.printf("The count of significant normal ASM Block (contains %d CpGs) is: %d\n",cpgNumInBlock, sigNormalCount);
		System.err.printf("The count of significant tumor ASM Block (contains %d CpGs) is: %d\n",cpgNumInBlock, sigTumorCount);
		System.err.printf("The count of overlap significant ASM block is: %d\n",sigOverlapCount);
		System.err.printf("The count of significant changed normal ASM Block (contains %d CpGs) is: %d\n",cpgNumInBlock, sigNormalChangeCount);
		System.err.printf("The count of significant changed tumor ASM Block (contains %d CpGs) is: %d\n",cpgNumInBlock, sigTumorChangeCount);
	}

}
