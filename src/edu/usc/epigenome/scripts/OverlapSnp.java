package edu.usc.epigenome.scripts;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.TreeSet;

import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;

public class OverlapSnp {

	/**
	 * @param args
	 */

	final private static String USAGE = "OverlapSnp chr snpfile1.txt snpfile2.txt ...";
	final private static String prefix = "methylCGsRich_ASM_";
	
	@Argument
	private List<String> stringArgs = new ArrayList<String>();
	
	public static void main(String[] args)
	throws Exception
	{
		new OverlapSnp().doMain(args);
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

		String chr = stringArgs.remove(0);
		String snpFileName1 = stringArgs.remove(0);
		String snpFileName2 = stringArgs.remove(0);
		BufferedReader br1 = new BufferedReader(new FileReader(snpFileName1));
		BufferedReader br2 = new BufferedReader(new FileReader(snpFileName2));
		String fn = prefix + chr + "_overlap_SNP_all_afterBaseQfilter.filterCNV" + ".txt";
		PrintWriter writer = new PrintWriter(new File(fn));
		String line;
		int count = 0;
		TreeSet<Integer> snpPosition1 = new TreeSet<Integer>();
		TreeSet<Integer> snpPosition2 = new TreeSet<Integer>();
		while( (line = br1.readLine()) != null){
			String[] tmpArray = line.split("\t");
			snpPosition1.add(Integer.parseInt(tmpArray[0]));
		}
		br1.close();
		while( (line = br2.readLine()) != null){
			String[] tmpArray = line.split("\t");
			snpPosition2.add(Integer.parseInt(tmpArray[0]));
		}
		br2.close();
		Iterator<Integer> it1 = snpPosition1.iterator();
		//Iterator<Integer> it2 = snpPosition2.iterator();
		while(it1.hasNext()){
			int tempSnp = it1.next();
			if(snpPosition2.contains(tempSnp)){
				writer.println(tempSnp);
				count++;
			}
		}
		writer.close();
		System.err.printf("The count of normal snp is: %d\n",snpPosition1.size());
		System.err.printf("The count of tumor snp is: %d\n",snpPosition2.size());
		System.err.printf("The count of overlap snp is: %d\n",count);
	}

}
