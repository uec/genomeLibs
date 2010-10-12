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
    protected String sample = "normal010310";
    //protected String sample = "test";
    
    //protected int minAlleleCount = 3;
    //protected double minAlleleFreq = 0.10;
    
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
			if (!chr.startsWith("chr")) chr = "chr" + chr;
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
		
		String sqlStatement = getSql(params);
		
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
			double pValue = fisherExact.getTwoTailedP(A_CReads,A_TReads,B_CReads,B_TReads);
			double logPValue = 0-Math.log10(pValue);
			String line = String.format("%d\t%d\t%f\t%f\t%d\t%d\t%d\t%d\t%c\t%c\t%c", methyCpg.alleleChromPos, methyCpg.chromPos, pValue, logPValue, A_CReads,A_TReads,B_CReads,B_TReads,methyCpg.getA_BaseUpperCase(), methyCpg.getB_BaseUpperCase(), methyCpg.getNextBaseRef());
			System.out.printf("%d\t%d\t%.2f\t%.2f\t%d\t%d\t%d\t%d\t%c\t%c\t%c\n", methyCpg.alleleChromPos, methyCpg.chromPos, pValue, logPValue, A_CReads,A_TReads,B_CReads,B_TReads,methyCpg.getA_BaseUpperCase(), methyCpg.getB_BaseUpperCase(), methyCpg.getNextBaseRef());
			outWriter.println(line);
			//System.err.printf("%d\t%d\t%d\t%d\t%d\t%d\t%c\t%c\t%c\n", methyCpg.alleleChromPos, methyCpg.chromPos, reads1[0], reads1[1],reads2[0],reads2[1],methyCpg.getA_BaseUpperCase(), methyCpg.getB_BaseUpperCase(), methyCpg.getNextBaseRef());

		}
		

		outWriter.close();
		
	}
	
	protected String getSql(MethylDbQuerier params) 
	throws Exception{
		String methTable = params.getMethylTable();
		//String methTable = params.methylTablePrefix;
		String sql = String.format("select * from %s WHERE ", methTable);
		sql += "ABaseRefUpperCase != '0'";
		sql += " AND BBaseRefUpperCase != '0'";
		sql += " AND (ACReads + ATReads >= 3)";
		sql += " AND (BCReads + BTReads >= 3)";
		sql += " AND (ACReads != 0 OR BCReads != 0)";
		sql += " AND (ATReads != 0 OR BTReads != 0)";
		sql += " AND (ACReads + ATReads)/totalReads >= 0.30";
		sql += " AND (BCReads + BTReads)/totalReads >= 0.30";
		if (Cpg)
			sql += " AND nextBaseRefUpperCase = 'G'";
		//sql += " GROUP BY chromPos "; // If you don't do this, you get multiple instances of the same CpG if it overlaps multiple features.
		sql += " ORDER BY chromPos,alleleChromPos ;";
		
		return sql;
	}
	
}


