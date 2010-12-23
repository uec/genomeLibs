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


public class MethylDbToPhase {

	/**
	 * @param args
	 */
	private static final String C_USAGE = "Use: MethylDbToPhase -tablePrefix " + MethylDbQuerier.DEFAULT_METHYL_TABLE_PREFIX + 
	" -windowSize 10 -GCH -sample s5 gff_file chr [startPos] [endPos]";
	//public static String connStr = null;
	protected static Connection cConn = null;
	public static String connStr = "jdbc:mysql://epifire2.epigenome.usc.edu/gnome_seq";
	protected int minOppoACount = 1;
    protected double minOppoAFreq = 0.10;
	//mysql_db_server: epifire2.epigenome.usc.edu
	
    @Option(name="-tablePrefix",usage="Prefix for DB table (default " + MethylDbQuerier.DEFAULT_METHYL_TABLE_PREFIX + ")")
    protected String tablePrefix = "methylCGsRich_gnome_seq_";
    @Option(name="-HCG",usage=" just withdarw HCpG sites")
    protected boolean Hcg = false;
    @Option(name="-GCH",usage=" just withdarw HCpG sites")
    protected boolean Gch = true;
    @Option(name="-sample",usage=" input the sample name: s5 or s6")
    protected String sample = "s6";
    @Option(name="-windowSize",usage=" input the window size: default is 10bp")
    protected int windowSize = 10;
    @Option(name="-rangeSize",usage=" input the range size: default is 2000bp")
    protected int rangeSize = 2000;
    
	// receives other command line parameters than options
	@Argument
	private List<String> arguments = new ArrayList<String>();
	
	public static void main(String[] args) 
	throws Exception{
		// TODO Auto-generated method stub
		new MethylDbToPhase().doMain(args);
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

			if(arguments.size() < 2 ) {
				System.err.println(C_USAGE);
				System.exit(1);
			}
			String gffFileName = arguments.get(0);
			BufferedReader br = new BufferedReader(new FileReader(gffFileName));
			String line;
			this.tablePrefix += sample;
			this.tablePrefix += "_";
			String fn1 = this.tablePrefix + "ctcf_table" + ".txt";
			PrintWriter outWriter = new PrintWriter(new File(fn1));
			String fn2 = this.tablePrefix + "ctcf_table_summary" + ".txt";
			PrintWriter outSumWriter = new PrintWriter(new File(fn2));
			Logger.getLogger(Logger.GLOBAL_LOGGER_NAME).setLevel(Level.INFO);
			if(connStr == null){
				BufferedReader hostBr = new BufferedReader(new FileReader("/home/uec-00/shared/production/database/mysql.host.txt"));
				String hostName = hostBr.readLine();
				connStr = "jdbc:mysql://" + hostName + "/gnome_seq";
			}
			setupDb(connStr);
			//int a = 0;
			while( (line = br.readLine()) != null){
				String[] tmpArray = line.split("\t");
				String chr = tmpArray[0];
				int chrStart = Integer.parseInt(tmpArray[3]); 
				int chrEnd = Integer.parseInt(tmpArray[4]);
				
				
				try {
					//step 4: create a statement
					String queryString = getSql(tablePrefix,chr,chrStart,chrEnd);
					Statement stmt = cConn.createStatement();
					ResultSet queryResult = stmt.executeQuery(queryString);
					double methyValue = queryResult.getDouble(0);
					queryResult.next();
					double readCoverage = queryResult.getDouble(0);
					outWriter.printf("%.2f\t%.2f\n", methyValue, readCoverage);
				}catch(SQLException ex) {
								System.err.println("Query: " + ex.getMessage());
				}
				
			}
			cleanupDb();
			outWriter.close();
			outSumWriter.close();

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
		String sql = String.format("select AVG(cReads/(cReads+tReads)), AVG(totalReads) from %s WHERE", methTable);
		sql += " AND totalReads != 0";
		sql += " AND (cReads+tReads) != 0";
		sql += " AND aReadsOpposite/totalReadsOpposite <= " + minOppoAFreq;
		sql += " AND aReadsOpposite <= " + minOppoACount;
		sql += " AND chromPos >= " + chrStart;
		sql += " AND chromPos <= " + chrEnd;
		if (Gch){
			sql += " AND nextBaseRefUpperCase != 'G'";
			sql += " AND preBaseRefUpperCase = 'G'";
		}
		else if(Hcg){
			sql += " AND nextBaseRefUpperCase = 'G'";
			sql += " AND preBaseRefUpperCase != 'G'";
		}
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
