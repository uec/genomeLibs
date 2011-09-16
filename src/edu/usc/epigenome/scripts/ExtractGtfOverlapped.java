package edu.usc.epigenome.scripts;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;

import org.biojava3.genome.parsers.gff.FeatureI;
import org.biojava3.genome.parsers.gff.FeatureList;
import org.biojava3.genome.parsers.gff.GFF3Reader;
import org.biojava3.genome.parsers.gff.GeneIDGFF2Reader;
import org.biojava3.genome.parsers.gff.Location;
import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;



public class ExtractGtfOverlapped {

	/**
	 * @param args
	 */
	
	
	final private static String USAGE = "ExtractGtfOverlapped [opt] Gtf1_whose_coordinate_will_output Gtf2_for_check GtfNameForOutput";
	
	@Option(name="-bothStrand",usage="check both strand overlap")
	protected boolean bothStrand = false;
	@Option(name="-excludeMode",usage="not overlap but export region that are not overlapped")
	protected boolean excludeMode = false;
	@Option(name="-overlapOrExcludeRegion",usage="the plus and minus region that are overlapped or excluded: -overlapOrExcludeRegion 2000 means +- 2kb excluded or overlapped")
	protected int overlapOrExcludeRegion = 0;
	
	@Argument
	private List<String> stringArgs = new ArrayList<String>();
	
	public static void main(String[] args) throws Exception {
		// TODO Auto-generated method stub
		//new ExtractGtfOverlapped().doMain(args);
		new ExtractGtfOverlapped(args);
	}

	
	public void doMain(String[] args)
	throws Exception {
		CmdLineParser parser = new CmdLineParser(this);
		parser.setUsageWidth(80);
		try
		{
			parser.parseArgument(args);
			if (stringArgs.size() < 3) throw new CmdLineException(USAGE);
			
			
		}
		catch (CmdLineException e)
		{
			System.err.println(e.getMessage());
			parser.printUsage(System.err);
			System.err.println();
			return;
		}


		
		
		String gtf1 = stringArgs.remove(0);
		String gtf2 = stringArgs.remove(0);
		String gtfOutput = stringArgs.remove(0);
		
		FeatureList gtf1FeatureList = GFF3Reader.read(gtf1);//my own made gtf should use GeneIDGFF2Reader; official gtf should use GFF3Reader
		FeatureList gtf2FeatureList = GFF3Reader.read(gtf2);
		
		Iterator<FeatureI> it1 = gtf1FeatureList.iterator();
		HashMap<String,FeatureList> gtf1Storage = new HashMap<String,FeatureList>();
		while(it1.hasNext()){
			FeatureI tmp = it1.next();
			FeatureList tmpList = new FeatureList();
			if(gtf1Storage.containsKey(tmp.seqname())){
				tmpList = gtf1Storage.get(tmp.seqname());
				tmpList.add(tmp);
			}
			else{
				tmpList.add(tmp);
			}
			gtf1Storage.put(tmp.seqname(), tmpList);
		}
		
		FeatureList overlappedGtfList = new FeatureList();
		Iterator<FeatureI> it2 = gtf2FeatureList.iterator();
		//HashMap<String,FeatureList> gtf2Storage = new HashMap<String,FeatureList>();
		while(it2.hasNext()){
			FeatureI tmp = it2.next();

			if(gtf1Storage.containsKey(tmp.seqname())){
				FeatureList tmpList = gtf1Storage.get(tmp.seqname());
				overlappedGtfList.add(tmpList.selectOverlapping(tmp.seqname(), tmp.location(), bothStrand));
			}


		}
		
		GeneIDGFF2Reader.write(overlappedGtfList, gtfOutput);
		
	}
	
	public ExtractGtfOverlapped(String[] args) throws Exception{
		CmdLineParser parser = new CmdLineParser(this);
		parser.setUsageWidth(80);
		try
		{
			parser.parseArgument(args);
			if (stringArgs.size() < 3) throw new CmdLineException(USAGE);
			
			
		}
		catch (CmdLineException e)
		{
			System.err.println(e.getMessage());
			parser.printUsage(System.err);
			System.err.println();
			return;
		}

		//String chr = stringArgs.remove(0);
		String genomicLocFileName = stringArgs.remove(0);
		String genomicLocRefFileName = stringArgs.remove(0);
		String fn = stringArgs.remove(0);
		BufferedReader genomicLocBr = new BufferedReader(new FileReader(genomicLocFileName));
		BufferedReader genomicLocRefBr = new BufferedReader(new FileReader(genomicLocRefFileName));
		
		PrintWriter outWriter = new PrintWriter(new File(fn));
		HashMap<String,List<Location>> genomicLocRefMap= new HashMap<String,List<Location>>();
	//	HashMap<String,List<Location>> genomicLocMap= new HashMap<String,List<Location>>();
		String line;
		int num = 0;
		while( (line = genomicLocRefBr.readLine()) != null ){
			if(num ==0){
				num++;
				continue;
			}
			String[] tmpArray = line.split("\t");
			if(genomicLocRefMap.containsKey(tmpArray[0])){
				List<Location> tempLoc = genomicLocRefMap.get(tmpArray[0]);
				Location thisLoc = new Location(Integer.parseInt(tmpArray[3]),Integer.parseInt(tmpArray[4]));
				tempLoc.add(thisLoc);
				//genomicLocMap.put(tmpArray[0], tempLoc);
				//String keyMark = tmpArray[0] + "~" + tmpArray[3] + "~" + tmpArray[4];
				genomicLocRefMap.put(tmpArray[0], tempLoc);
			}
			else{
				List<Location> tempLoc = new ArrayList<Location>();
				Location thisLoc = new Location(Integer.parseInt(tmpArray[3]),Integer.parseInt(tmpArray[4]));
				tempLoc.add(thisLoc);
			//	genomicLocMap.put(tmpArray[0], tempLoc);
				//String keyMark = tmpArray[0] + "~" + tmpArray[3] + "~" + tmpArray[4];
				genomicLocRefMap.put(tmpArray[0], tempLoc);
			}
		}
		genomicLocRefBr.close();
		num=0;
		while( (line = genomicLocBr.readLine()) != null ){
			if(num ==0){
				num++;
				continue;
			}
			String[] tmpArray = line.split("\t");
			int start = Integer.parseInt(tmpArray[3]) - overlapOrExcludeRegion;
			int end = Integer.parseInt(tmpArray[4]) + overlapOrExcludeRegion; //do not check overflow contig length
			if(start <= 0){
				start = Integer.parseInt(tmpArray[3]);
			}
			Location checkLoc = new Location(start,end);
			if(genomicLocRefMap.containsKey(tmpArray[0])){
				
				List<Location> thisLoc = genomicLocRefMap.get(tmpArray[0]);
				Iterator<Location> It = thisLoc.iterator();
				boolean overlapped = false;
				while(It.hasNext()){
					Location tempLoc = It.next();
					if(excludeMode){
						if(tempLoc.overlaps(checkLoc)){
							//	String keySearch = tmpArray[0] + "~" + tempLoc.getBegin() + "~" + tempLoc.getEnd();
								//String content = genomicLocRefMap.get(keySearch);
								overlapped = true;
								//System.out.println(tempLoc);
								//System.out.println(checkLoc.toString());
								//System.out.println(tempLoc.contains(checkLoc));
								break;
						}
					}
					else{
						if(tempLoc.overlaps(checkLoc)){
							//	String keySearch = tmpArray[0] + "~" + tempLoc.getBegin() + "~" + tempLoc.getEnd();
								//String content = genomicLocRefMap.get(keySearch);
								outWriter.println(line);
								//System.out.println(tempLoc);
								//System.out.println(checkLoc.toString());
								//System.out.println(tempLoc.contains(checkLoc));
								break;
						}
					}
				}
				if(excludeMode){
					if(!overlapped){
						outWriter.println(line);
					}
				}
			}
		}
		outWriter.close();
		genomicLocBr.close();
	}
		
}
