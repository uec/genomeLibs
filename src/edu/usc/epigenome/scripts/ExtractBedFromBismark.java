package edu.usc.epigenome.scripts;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;

import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;


//notice, bed file are 0 based coordinate!!! The genome browser region "chr1:1-1000" would be described in a BED record as "chr1 0 1000" with the start coordinate 
//being one smaller and the end coordinate being the same, describing the half-closed half-open interval [0,1000) of length 1000bp starting at base 0.
public class ExtractBedFromBismark {

	private static final String C_USAGE = "Use: ExtractBedFromBismark inputBismarkFile";
	
	// receives other command line parameters than options
	@Argument
	private List<String> arguments = new ArrayList<String>();
	
	/**
	 * @param args
	 */
	
	@Option(name="-coverage",usage=" minimum CT reads coverage")
    protected int coverage = 1;
	@Option(name="-minMethyPercent",usage=" minimum Methy Percent")
    protected int minMethyPercent = 0;
	
	
	
	public static void main(String[] args) 
	throws Exception{
		// TODO Auto-generated method stub
		new ExtractBedFromBismark().doMain(args);
	}

	public void doMain(String[] args)
	throws Exception
	{
		CmdLineParser parser = new CmdLineParser(this);
		// if you have a wider console, you could increase the value;
		// here 80 is also the default
		
		parser.setUsageWidth(80);
		try
		{
			parser.parseArgument(args);

			if(arguments.size() < 1 ) {
				System.err.println(C_USAGE);
				System.exit(1);
			}
			String bisFileName = arguments.get(0);
			BufferedReader br = new BufferedReader(new FileReader(bisFileName));
			String line;
			String fn = bisFileName + ".bed";
			PrintWriter outWriter = new PrintWriter(new File(fn));
			HashMap<String,Integer[]> cytosine = new HashMap<String,Integer[]>();

			while( (line = br.readLine()) != null){
				if(line.startsWith("Bismark"))
					continue;
				String[] tmpArray = line.split("\t");
				String position = tmpArray[2] + "~" + tmpArray[3]; 
				String methyState = tmpArray[1];
				String methyState2 = tmpArray[4];
				Integer[] methyValue = new Integer[2]; // [0] = methylated, [1] = total
				if(cytosine.containsKey(position)){
					methyValue = cytosine.get(position);
				}
				else{
					methyValue[0] = 0;
					methyValue[1] = 0;
				}
				
				boolean passValidate = validateMethylationCall(methyState,methyState2);
				if(passValidate){
					if(methyState.equalsIgnoreCase("+")){
						methyValue[0]++;
						methyValue[1]++;
					}
					else{
						
						methyValue[1]++;
					}
					cytosine.put(position, methyValue);
				}
				
			}
			br.close();
			
			Iterator<String> It = cytosine.keySet().iterator();
			while(It.hasNext()){
				String key = It.next();
				String[] tmpArray = key.split("~");
				String chr = tmpArray[0];
				int start = Integer.parseInt(tmpArray[1]);//bed file are 0 based coordinate
				int end = start + 1;
				Integer[] tempMethyValue = cytosine.get(key);
				double methyPercentage = 0;
				if(tempMethyValue[1] >= coverage){
					methyPercentage = (double)tempMethyValue[0]/(double)tempMethyValue[1];
					if(methyPercentage >= minMethyPercent){
						outWriter.println(chr + "\t" + start + "\t" + end + "\t" + methyPercentage);
					}
				}
			}
			outWriter.close();
			
		}
		catch (CmdLineException e)
		{
			System.err.println(e.getMessage());
			System.err.println(C_USAGE);
			// print the list of available options
			parser.printUsage(System.err);
			System.err.println();
			return;
		}	
	}
	
	private boolean validateMethylationCall(String methyState, String methyState2){
		boolean passValidate = false;
		if(methyState2.equalsIgnoreCase("z")){
			passValidate = validateCpgMethylationCall(methyState,methyState2);
		}
		else{
			passValidate = validateNonCpgMethylationCall(methyState,methyState2);
		}
		
		return passValidate;
	}
	
	private boolean validateCpgMethylationCall(String methyState, String methyState2){
		boolean passValidate = false;
		if(methyState2.equalsIgnoreCase("z") && methyState.equalsIgnoreCase("-")){
			passValidate = true;
			return passValidate;
		}
		if(methyState2.equalsIgnoreCase("Z") && methyState.equalsIgnoreCase("+")){
			passValidate = true;
			return passValidate;
		}
		
		return passValidate;
	}
	
	private boolean validateNonCpgMethylationCall(String methyState, String methyState2){
		boolean passValidate = false;
		if(methyState2.equalsIgnoreCase("c") && methyState.equalsIgnoreCase("-")){
			passValidate = true;
			return passValidate;
		}
		if(methyState2.equalsIgnoreCase("C") && methyState.equalsIgnoreCase("+")){
			passValidate = true;
			return passValidate;
		}
		if(methyState2.equalsIgnoreCase("x") && methyState.equalsIgnoreCase("-")){
			passValidate = true;
			return passValidate;
		}
		if(methyState2.equalsIgnoreCase("X") && methyState.equalsIgnoreCase("+")){
			passValidate = true;
			return passValidate;
		}
		if(methyState2.equalsIgnoreCase("h") && methyState.equalsIgnoreCase("-")){
			passValidate = true;
			return passValidate;
		}
		if(methyState2.equalsIgnoreCase("H") && methyState.equalsIgnoreCase("+")){
			passValidate = true;
			return passValidate;
		}

		return passValidate;
	}
}
