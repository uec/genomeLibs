package edu.usc.epigenome.scripts;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.PrintWriter;
import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.TreeSet;
import java.util.TreeMap;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.lang.Number;

import net.sf.samtools.SAMFileReader;

import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;

import com.sun.org.apache.xalan.internal.xsltc.compiler.Pattern;

import edu.usc.epigenome.genomeLibs.FisherExactTest;
import edu.usc.epigenome.genomeLibs.GoldAssembly;

import edu.usc.epigenome.genomeLibs.MethylDb.Cytosine;
import edu.usc.epigenome.genomeLibs.MethylDb.CytosineIterator;
import edu.usc.epigenome.genomeLibs.MethylDb.MethylDbQuerier;


public class MethylDbToPhaseExtend {

	/**
	 * @param args
	 */
	private static final String C_USAGE = "Use: MethylDbToPhase -tablePrefix " + "methylCGsRich_gnome_seq_merge_s5_s6_" + "-sample ctcf" +
	" -windowSize 20 -GCH gff_file";
	//public static String connStr = null;
	protected static Connection cConn = null;
	public static String connStr = "jdbc:mysql://epifire2.epigenome.usc.edu/gnome_seq";
	protected int minOppoACount = 1;
    protected double minOppoAFreq = 0.10;
	//mysql_db_server: epifire2.epigenome.usc.edu
	
    @Option(name="-tablePrefix",usage="Prefix for DB table (default " + MethylDbQuerier.DEFAULT_METHYL_TABLE_PREFIX + ")")
    protected String tablePrefix = "methylCGsRich_gnome_seq_merge_s5_s6_";
    @Option(name="-WCG",usage=" just withdarw A/TCpG sites")
    protected boolean Wcg = false;
    @Option(name="-CCG",usage=" just withdarw CCpG sites")
    protected boolean Ccg = false;
    @Option(name="-GCH",usage=" just withdarw HCpG sites")
    protected boolean Gch = false;
    @Option(name="-sample",usage=" input the sample name: ctcf")
    protected String sample = "ctcf";
    @Option(name="-windowSize",usage=" input the window size: default is 20bp")
    protected int windowSize = 20;
    @Option(name="-rangeSize",usage=" input the range size: default is 2000bp")
    protected int rangeSize = 2000;
    @Option(name="-alignmentType",usage="alignment type(1: 5' end; 2: 3' end; 3: center) the default is 3")
    protected int alignmentType = 3;
    @Option(name="-hpcc",usage=" use hpcc database")
    protected boolean hpcc = false;
    
	// receives other command line parameters than options
	@Argument
	private List<String> arguments = new ArrayList<String>();
	
	public static void main(String[] args) 
	throws Exception{
		// TODO Auto-generated method stub
		new MethylDbToPhaseExtend().doMain(args);
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
			BufferedReader br = new BufferedReader(new FileReader(gffFileName));
			String line;
			
			String fn1 = this.tablePrefix + sample;
			String fn2 = this.tablePrefix + sample;
			if(Gch){
				fn1 += "_gch_table.txt";
				fn2 += "_gch_table.addtional.txt";
			}
			else if(Wcg){
				fn1 += "_wcg_table.txt";
				fn2 += "_wcg_table.addtional.txt";
			}
			else if(Ccg){
				fn1 += "_ccg_table.txt";
				fn2 += "_ccg_table.addtional.txt";
			}
			else{
				fn1 += "_table.txt";
				fn2 += "_table.addtional.txt";
			}
			PrintWriter outWriter = new PrintWriter(new File(fn1));
			//additional file to store c reads, t reads number detail
			PrintWriter outWriter2 = new PrintWriter(new File(fn2));

			Logger.getLogger(Logger.GLOBAL_LOGGER_NAME).setLevel(Level.INFO);
			
			if(connStr == null || hpcc){
				BufferedReader hostBr = new BufferedReader(new FileReader("/home/uec-00/shared/production/database/mysql.host.txt"));
				String hostName = hostBr.readLine();
				connStr = "jdbc:mysql://" + hostName + "/gnome_seq";
			}
			setupDb(connStr);

			int recCounter = 0;
			while( (line = br.readLine()) != null){
				if(!line.startsWith("chr"))
					continue;
				String[] tmpArray = line.split("\t");
				String chr = tmpArray[0];
				int chrStart = Integer.parseInt(tmpArray[3]); 
				int chrEnd = Integer.parseInt(tmpArray[4]);
				boolean negStrand = tmpArray[6].charAt(0) == '-' ? true : false;
				int motifCenter = 0;
				if(negStrand){
					if(alignmentType == 3){
						motifCenter = (chrEnd + chrStart)/2;
					}
					else if(alignmentType == 2){
						motifCenter = chrStart;
					}
					else if(alignmentType ==1){
						motifCenter = chrEnd;
					}
					else{
						System.err.println("alignment type error!");
					}
				}
				else{
					if(alignmentType == 3){
						motifCenter = (chrEnd + chrStart)/2;
					}
					else if(alignmentType == 2){
						motifCenter = chrEnd;
					}
					else if(alignmentType ==1){
						motifCenter = chrStart;
					}
					else{
						System.err.println("alignment type error!");
					}
				}
				
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
				System.err.println("motifCenter: " + motifCenter + " negStrand: " + negStrand);
				recCounter++;
				if(negStrand){
					for(int windowStart = rangeEnd - windowSize, windowEnd = rangeEnd, i = 0; windowStart >= rangeStart; windowStart -= windowSize, windowEnd -= windowSize, i++){
						try {
							//step 4: create a statement
							String queryString = getSql(tablePrefix,chr,windowStart,windowEnd);
							Statement stmt = cConn.createStatement();
							ResultSet queryResult = stmt.executeQuery(queryString);
							queryResult.next();
							
							double methyValue = queryResult.getDouble(1);
							int numC = queryResult.getInt(2);
							double cReads = queryResult.getDouble(3);
							double tReads = queryResult.getDouble(4);

							if(numC != 0){

								outWriter.printf("%.2f\t", methyValue);
								outWriter2.printf("%d\t%.2f\t%.2f\t", numC, cReads, tReads);

							}
							else{
								outWriter.printf("NA\t");
								outWriter2.printf("NA\tNA\tNA\t");
							}
							queryResult.close();
							
						}catch(SQLException ex) {
										System.err.println("Query: " + ex.getMessage());
						}
					}
				}
				else{
					for(int windowStart = rangeStart, windowEnd = rangeStart + windowSize, i = 0; windowEnd <= rangeEnd; windowStart+=windowSize, windowEnd+=windowSize, i++){
						try {
							//step 4: create a statement
							String queryString = getSql(tablePrefix,chr,windowStart,windowEnd);
							Statement stmt = cConn.createStatement();
							ResultSet queryResult = stmt.executeQuery(queryString);
							queryResult.next();
							
							double methyValue = queryResult.getDouble(1);
							
							
							int numC = queryResult.getInt(2);
							double cReads = queryResult.getDouble(3);
							double tReads = queryResult.getDouble(4);

							if(numC != 0){

								outWriter.printf("%.2f\t", methyValue);
								outWriter2.printf("%d\t%.2f\t%.2f\t", numC, cReads, tReads);

							}
							else{
								outWriter.printf("NA\t");
								outWriter2.printf("NA\tNA\tNA\t");
							}
							queryResult.close();
							
						}catch(SQLException ex) {
										System.err.println("Query: " + ex.getMessage());
						}
					}
				}
				
				outWriter.printf("\n");
				outWriter2.printf("\n");

				
			}
			cleanupDb();
			outWriter.close();
			outWriter2.close();		

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
	
	protected String getSql(String tablePrefix, String chr, int chrStart, int chrEnd) 
	throws Exception{
		String methTable = tablePrefix + chr;
		String sql = String.format("select AVG(cReads/(cReads+tReads)), COUNT(*), AVG(cReads), AVG(tReads) from %s WHERE", methTable);
		sql += " totalReads != 0";
		sql += " AND (cReads+tReads) != 0";
		sql += " AND aReadsOpposite/totalReadsOpposite <= " + minOppoAFreq;
		sql += " AND aReadsOpposite <= " + minOppoACount;
		sql += " AND chromPos >= " + chrStart;
		sql += " AND chromPos <= " + chrEnd;
		if (Gch){
			sql += " AND nextBaseRefUpperCase != 'G'";
			sql += " AND preBaseRefUpperCase = 'G'";
		}
		else if(Wcg){
			sql += " AND nextBaseRefUpperCase = 'G'";
			sql += " AND (preBaseRefUpperCase = 'A' OR preBaseRefUpperCase = 'T')";
		}
		else if(Ccg){
			sql += " AND nextBaseRefUpperCase = 'G'";
			sql += " AND preBaseRefUpperCase = 'C'";
		}
		//System.err.println(sql);
		return sql;
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
