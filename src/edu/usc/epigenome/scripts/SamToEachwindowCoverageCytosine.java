package edu.usc.epigenome.scripts;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.PrintWriter;
import java.sql.Connection;
import java.util.ArrayList;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;

import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.util.CloseableIterator;

import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;

import edu.usc.epigenome.genomeLibs.MethylDb.MethylDbQuerier;

public class SamToEachwindowCoverageCytosine {
	/**
	 * @param args
	 */
	private static final String C_USAGE = "Use: SamToEachwindowCoverageCytosine -tablePrefix " + "methylCGsRich_gnome_seq_merge_s5_s6_" + "-sample ctcf" +
	" -windowSize 10 -minMapQ 0 gff_file bam_file";
	//public static String connStr = null;
	//protected static Connection cConn = null;
	//public static String connStr = "jdbc:mysql://epifire2.epigenome.usc.edu/gnome_seq";
	//protected int minOppoACount = 2;
    //protected double minOppoAFreq = 0.20;
	//mysql_db_server: epifire2.epigenome.usc.edu
	
    @Option(name="-tablePrefix",usage="Prefix for DB table (default " + MethylDbQuerier.DEFAULT_METHYL_TABLE_PREFIX + ")")
    protected String tablePrefix = "methylCGsRich_window_coverage_";
    //@Option(name="-HCG",usage=" just withdarw HCpG sites")
    //protected boolean Hcg = false;
    //@Option(name="-GCH",usage=" just withdarw GpCH sites")
    //protected boolean Gch = false;
    @Option(name="-sample",usage=" input the sample name: ctcf")
    protected String sample = "ctcf";
    @Option(name="-windowSize",usage=" input the window size: default is 10bp")
    protected int windowSize = 20;
    @Option(name="-rangeSize",usage=" input the range size: default is 2000bp")
    protected int rangeSize = 2000;
    @Option(name="-minMapQ",usage="minimum mapping quality (default 30)")
	protected int minMapQ = 0;
    
	// receives other command line parameters than options
	@Argument
	private List<String> arguments = new ArrayList<String>();
	
	public static void main(String[] args) 
	throws Exception{
		// TODO Auto-generated method stub
		new SamToEachwindowCoverageCytosine().doMain(args);
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
			String gffFileName = arguments.get(0);
			String fn1 = arguments.get(1);
			BufferedReader br = new BufferedReader(new FileReader(gffFileName));
			String line;
			//this.tablePrefix += sample;
			//this.tablePrefix += "_";
			File inputSamOrBamFile = new File(fn1);
			final SAMFileReader inputSam = new SAMFileReader(inputSamOrBamFile);
			inputSam.setValidationStringency(SAMFileReader.ValidationStringency.SILENT);
			
			String fn2 = this.tablePrefix + sample + "_table_coverage" + ".txt";
			PrintWriter outCovWriter = new PrintWriter(new File(fn2));
			//String fn2 = this.tablePrefix + sample + "_table_summary" + ".txt";
			//PrintWriter outSumWriter = new PrintWriter(new File(fn2));
			Logger.getLogger(Logger.GLOBAL_LOGGER_NAME).setLevel(Level.INFO);
			//int a = 0;
			//int gffLineNumber = 0;
			//TreeMap<Integer,List<Double>> methySum = new TreeMap<Integer,List<Double>>();
			int recCounter = 0;
			while( (line = br.readLine()) != null){
				String[] tmpArray = line.split("\t");
				String chr = tmpArray[0];
				int chrStart = Integer.parseInt(tmpArray[3]); 
				int chrEnd = Integer.parseInt(tmpArray[4]);
				boolean negStrand = tmpArray[6].charAt(0) == '-' ? true : false;
				int motifCenter = (chrEnd + chrStart)/2;
				int rangeStart = motifCenter - rangeSize;
				int rangeEnd = motifCenter + rangeSize;
				if(rangeStart < 0 || rangeEnd < 0 || rangeStart >= rangeEnd){
					System.err.println("motifCenter < 0");
				}
				if ((recCounter % 100)==0)
				{
					System.err.printf("On new record #%d\n",recCounter);
					System.gc();
				}
				recCounter++;
				
				if(negStrand){
					for(int windowStart = rangeEnd - windowSize, windowEnd = rangeEnd, i = 0; windowStart >= rangeStart; windowStart -= windowSize, windowEnd -= windowSize, i++){
						int coverage = 0;
						CloseableIterator<SAMRecord> chrIt = inputSam.queryOverlapping(chr, windowStart, windowEnd);
						record:while (chrIt.hasNext())
						{
							SAMRecord samRecord = chrIt.next();
							int mapQual = samRecord.getMappingQuality();
							boolean unmapped = samRecord.getReadUnmappedFlag();
							if (unmapped || (mapQual <= minMapQ))
							{
								continue record;
							}
							coverage++;
						}
						outCovWriter.printf("%d\t", coverage);
						chrIt.close();
					}
				}
				else{
					for(int windowStart = rangeStart, windowEnd = rangeStart + windowSize, i = 0; windowEnd <= rangeEnd; windowStart+=windowSize, windowEnd+=windowSize, i++){
						int coverage = 0;
						CloseableIterator<SAMRecord> chrIt = inputSam.queryOverlapping(chr, windowStart, windowEnd);
						record:while (chrIt.hasNext())
						{
							SAMRecord samRecord = chrIt.next();
							int mapQual = samRecord.getMappingQuality();
							boolean unmapped = samRecord.getReadUnmappedFlag();
							if (unmapped || (mapQual <= minMapQ))
							{
								continue record;
							}
							coverage++;
						}
						outCovWriter.printf("%d\t", coverage);
						chrIt.close();
					}
				}
				outCovWriter.printf("\n");
				
			}
			outCovWriter.close();
			inputSam.close();
			/*Iterator<List<Double>> windowIt = methySum.values().iterator();
			double methyValueSum = 0;
			int count = 0;
			while(windowIt.hasNext()){
				List<Double> windowValue = windowIt.next();
				Iterator<Double> methyIt = windowValue.iterator();
				count = windowValue.size();
				while(methyIt.hasNext()){
					methyValueSum += methyIt.next();
				}
				outSumWriter.printf("%.2f\t",methyValueSum/(double)count);
			}*/
			//outSumWriter.printf("\n");
			//outSumWriter.close();
			

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
	
}
