package edu.usc.epigenome.scripts;

import java.io.File;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.TreeSet;
import java.util.TreeMap;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.lang.Number;

import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;

import com.sun.org.apache.xalan.internal.xsltc.compiler.Pattern;

import edu.usc.epigenome.genomeLibs.GoldAssembly;

import edu.usc.epigenome.genomeLibs.MethylDb.Cytosine;
import edu.usc.epigenome.genomeLibs.MethylDb.CytosineIterator;
import edu.usc.epigenome.genomeLibs.MethylDb.MethylDbQuerier;
import edu.usc.epigenome.genomeLibs.FisherExactTest;;



public class MethylDbToFisherExactLogpASM {

	/**
	 * @param args
	 */
	private static final String C_USAGE = "Use: MethylDbToChisquareLogpASM -tablePrefix " + MethylDbQuerier.DEFAULT_METHYL_TABLE_PREFIX + 
	" CpG sample chr [startPos] [endPos]";
	public static String connStr = "jdbc:mysql://hpc2721/asm_cr";
	//public static String connStr = "jdbc:mysql://epifire2.epigenome.usc.edu/asm_cr";
	//mysql_db_server: epifire2.epigenome.usc.edu
	
    @Option(name="-tablePrefix",usage="Prefix for DB table (default " + MethylDbQuerier.DEFAULT_METHYL_TABLE_PREFIX + ")")
    protected String tablePrefix = "methylCGsRich_ASM_";
    //protected String tablePrefix = "methylCGsRich_test_";
    @Option(name="-CpG",usage=" just withdarw CpG sites or all of the cytosine sites")
    protected boolean Cpg = true;
    @Option(name="-sample",usage=" input the sample name: normal010310 or tumor011010")
    //protected String sample = "normal010310";
    protected String sample = "test";
    
    protected int minAlleleCount = 3;
    protected double minAlleleFreq = 0.30;
    protected int minOppoACount = 1;
    protected double minOppoAFreq = 0.10;
    
	// receives other command line parameters than options
	@Argument
	private List<String> arguments = new ArrayList<String>();
	
	public static void main(String[] args) 
	throws Exception{
		// TODO Auto-generated method stub
		new MethylDbToFisherExactLogpASM().doMain(args);
	}

	public void doMain(String[] args)
	throws Exception
	{
		CmdLineParser parser = new CmdLineParser(this);
		// if you have a wider console, you could increase the value;
		// here 80 is also the default
		int chrStart = -1, chrEnd = -1;
		String chr;
		parser.setUsageWidth(80);
		try
		{
			parser.parseArgument(args);

			if(arguments.size() < 1 ) {
				System.err.println(C_USAGE);
				System.exit(1);
			}

			chr = arguments.get(0);
			if (!chr.startsWith("chr")) {
				chr = chr.toUpperCase();
				chr = "chr" + chr;
			}
			System.err.println(chr);
			if (arguments.size() > 1)
			{
				chrStart = Integer.parseInt(arguments.get(1));
				chrEnd = Integer.parseInt(arguments.get(2));
			}

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
		this.tablePrefix += sample;
		this.tablePrefix += "_";
		
		Logger.getLogger(Logger.GLOBAL_LOGGER_NAME).setLevel(Level.INFO);
		
		MethylDbQuerier params = new MethylDbQuerier();
		if (this.tablePrefix != null) params.methylTablePrefix = this.tablePrefix;
		params.setMinCTreads(0);
		params.setUseNonconversionFilter(false);

		if (chrStart<0)
		{
			chrStart = 1;
			chrEnd = GoldAssembly.chromLengthStatic(chr, "hg18");
			//System.out.println(chrEnd);
		}
		
		params.clearRangeFilters();
		params.addRangeFilter(chr, chrStart, chrEnd);
		
		String fn = this.tablePrefix + chr + "_pValue" + ".txt";
		PrintWriter outWriter = new PrintWriter(new File(fn));
		
		String sqlStatement = getSql(params, chr);
		
		CytosineIterator it = new CytosineIterator(params,connStr,sqlStatement);
		TreeMap<Integer,Cytosine> totalReadsSum = new TreeMap<Integer,Cytosine>();
		//get rid of duplicate cytosine position. 
		while(it.hasNext()){
			Cytosine methyCpg = it.next(true);
			
			if(totalReadsSum.containsKey(methyCpg.chromPos)){
				if(methyCpg.totalReads > totalReadsSum.get(methyCpg.chromPos).totalReads){
					totalReadsSum.put(methyCpg.chromPos, methyCpg);
				}
				else{
					continue;
				}
			}
			else{
				totalReadsSum.put(methyCpg.chromPos, methyCpg);
			}
			
						
		}
		
		Iterator<Cytosine> cytocineOutput = totalReadsSum.values().iterator();
		while (cytocineOutput.hasNext())
		{
			Cytosine methyCpg = cytocineOutput.next();
			FisherExactTest fisherExact = new FisherExactTest(100);
			int A_CReads = methyCpg.getA_CReads();
			int A_TReads = methyCpg.getA_TReads();
			int B_CReads = methyCpg.getB_CReads();
			int B_TReads = methyCpg.getB_TReads();
			if(methyCpg.totalReadsOpposite <= 10 && methyCpg.aReadsOpposite > minOppoACount )
				 continue;
			double[] logExpect = logExpectation(A_CReads,A_TReads,B_CReads,B_TReads);
			double pValue = fisherExact.getTwoTailedP(A_CReads,A_TReads,B_CReads,B_TReads);
			double logPValue = 0-Math.log(pValue)/Math.log(2);
			 
			String line = String.format("%d\t%d\t%f\t%f\t%f\t%f\t%d\t%d\t%d\t%d\t%c\t%c\t%c", methyCpg.alleleChromPos, methyCpg.chromPos, pValue, logPValue, logExpect[0], logExpect[1], A_CReads,A_TReads,B_CReads,B_TReads,methyCpg.getA_BaseUpperCase(), methyCpg.getB_BaseUpperCase(), methyCpg.getNextBaseRef());
			//System.out.printf("%d\t%d\t%.2f\t%.2f\t%d\t%d\t%d\t%d\t%c\t%c\t%c\n", methyCpg.alleleChromPos, methyCpg.chromPos, pValue, logPValue, A_CReads,A_TReads,B_CReads,B_TReads,methyCpg.getA_BaseUpperCase(), methyCpg.getB_BaseUpperCase(), methyCpg.getNextBaseRef());
			outWriter.println(line);
			//System.err.printf("%d\t%d\t%d\t%d\t%d\t%d\t%c\t%c\t%c\n", methyCpg.alleleChromPos, methyCpg.chromPos, reads1[0], reads1[1],reads2[0],reads2[1],methyCpg.getA_BaseUpperCase(), methyCpg.getB_BaseUpperCase(), methyCpg.getNextBaseRef());

		}
		

		outWriter.close();
		
	}
	
	protected String getSql(MethylDbQuerier params, String chr) 
	throws Exception{
		//String methTable = params.getMethylTable();
		String methTable = tablePrefix + chr;
		//String methTable = params.methylTablePrefix;
		String sql = String.format("select * from %s WHERE ", methTable);
		sql += "ABaseRefUpperCase != '0'";
		sql += " AND BBaseRefUpperCase != '0'";
		
		sql += " AND (ACReads != 0 OR BCReads != 0)";
		sql += " AND (ATReads != 0 OR BTReads != 0)";
		//sql += " AND (ACReads + ATReads)/totalReads >= " + minAlleleFreq;
		//sql += " AND (BCReads + BTReads)/totalReads >= " + minAlleleFreq;
		sql += " AND (ACReads + ATReads >= " + minAlleleCount + ")";
		sql += " AND (BCReads + BTReads >= " + minAlleleCount + ")";
		sql += " AND aReadsOpposite/totalReadsOpposite <= " + minOppoAFreq;
		if (Cpg){
			sql += " AND nextBaseRefUpperCase = 'G'";
		}
			
		//sql += " GROUP BY chromPos "; // If you don't do this, you get multiple instances of the same CpG if it overlaps multiple features.
		sql += " ORDER BY chromPos,alleleChromPos ;";
		
		return sql;
	}
	
	protected double[] logExpectation(double A_CReads, double A_TReads, double B_CReads, double B_TReads){
		double biggerMethy = A_CReads > B_CReads ? A_CReads + A_TReads : B_CReads + B_TReads;
		double biggerUnmethy = A_TReads > B_TReads ? A_CReads + A_TReads : B_CReads + B_TReads;
		double methyExpect = (biggerMethy/(A_CReads + B_CReads + A_TReads + B_TReads)) * ((A_CReads + B_CReads)/(A_CReads + B_CReads + A_TReads + B_TReads)) * (A_CReads + B_CReads + A_TReads + B_TReads);
		double unmethyExpect = (biggerUnmethy/(A_CReads + B_CReads + A_TReads + B_TReads)) * ((A_TReads + B_TReads)/(A_CReads + B_CReads + A_TReads + B_TReads)) * (A_CReads + B_CReads + A_TReads + B_TReads);
		double[] logExpect = new double[2];
		logExpect[0] = (A_CReads > B_CReads ? A_CReads : B_CReads)/methyExpect;
		logExpect[1] = (A_TReads > B_TReads ? A_TReads : B_TReads)/unmethyExpect;
		logExpect[0] = Math.log(logExpect[0])/Math.log(2);
		logExpect[1] = Math.log(logExpect[1])/Math.log(2);
		return logExpect;
	}
	
}


