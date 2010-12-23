package edu.usc.epigenome.genomeLibs.MethylDb;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.SortedMap;
import java.util.TreeMap;
import java.sql.*;

import org.biojava.bio.seq.StrandedFeature;

import edu.usc.epigenome.genomeLibs.GoldAssembly;
import edu.usc.epigenome.genomeLibs.MethylDb.CytosineIterator;
import edu.usc.epigenome.genomeLibs.MethylDb.MethylDbQuerier;

public class Cytosine implements Comparable, Cloneable {

	private static final double DEFAULT_FAIL_MIN_FRAC_A = 0.2;
	private static final int DEFAULT_FAIL_MIN_ABS_A = 2;
	
	// Object vars
	public int chromPos = 0;
	public boolean negStrand = false;
	public short totalReads = 0;
	public short cReads = 0;
	public short cReadsNonconversionFilt = 0;
	public short tReads = 0;
	public short agReads = 0;
	public short totalReadsOpposite = 0;
	public short aReadsOpposite = 0;
	
	public short nextBaseTotalReads = 0;
	public short nextBaseGreads = 0;
	
	public short preBaseTotalReads = 0;
	public short preBaseGreads = 0;
	
	public int alleleChromPos = -1;
	protected char A_BaseUpperCase = '0';
	protected char B_BaseUpperCase = '0';
	public short A_CReads = 0; 
	public short B_CReads = 0;
	public short A_TReads = 0;
	public short B_TReads = 0;
	
	// This is actually for non-cpgs
	protected char nextBaseRefUpperCase = '0';
	protected char preBaseRefUpperCase = '0';
	
	protected double methyDens = Double.NaN;
	protected double methyA_Dens = Double.NaN;
	protected double methyB_Dens = Double.NaN;
	
	// This is for weighted averages
	protected double gchWeight = -1.00;
	protected double gcgWeight = -1.00;
	protected double hcgWeight = -1.00;

	protected static Connection conn = null;

	

	/**
	 * 
	 */
	public Cytosine() {
		super();
	}

	/**
	 * @param chromPos
	 * @param negStrand
	 * @param totalReads
	 * @param cReads
	 * @param cReadsNonconversionFilt
	 * @param tReads
	 * @param agReads
	 * @param totalReadsOpposite
	 * @param aReadsOpposite
	 */
	public Cytosine(int chromPos, boolean negStrand) {
		super();
		this.chromPos = chromPos;
		this.negStrand = negStrand;
	}
	
	public Cytosine(int chromPos, boolean negStrand, int alleleChromPos) {
		super();
		this.chromPos = chromPos;
		this.negStrand = negStrand;
		this.alleleChromPos = alleleChromPos;
	}
	
	public Cytosine(int chromPos, boolean negStrand, short totalReads, short cReads,
			short cReadsNonconversionFilt, short tReads, short agReads,
			short totalReadsOpposite, short aReadsOpposite, String preBaseRefUpperCase, String nextBaseRefUpperCase, double fracMeth, int gchWeight,  int gcgWeight, int hcgWeight ) {
		super();
		this.chromPos = chromPos;
		this.negStrand = negStrand;
		this.totalReads = totalReads;
		this.cReads = cReads;
		this.cReadsNonconversionFilt = cReadsNonconversionFilt;
		this.tReads = tReads;
		this.agReads = agReads;
		this.totalReadsOpposite = totalReadsOpposite;
		this.aReadsOpposite = aReadsOpposite;
		
		setPreBaseRef(preBaseRefUpperCase);
		setNextBaseRef(nextBaseRefUpperCase);
		
		this.methyDens = fracMeth;
		
		// CpG weight has some issues.  For instance at large gaps you get huge values.  Remove these
		if (gchWeight < 0)
		{
			System.err.printf("Got negative gchWeight: pos=%d, weight=%d\n",chromPos,gchWeight);
			gchWeight = 100;
		}
		else if (gchWeight > 5000)
		{
			System.err.printf("Got extra-large gchWeight: pos=%d, weight=%d\n",chromPos,gchWeight);
			gchWeight = 100;
		}
		
		if (gcgWeight < 0)
		{
			System.err.printf("Got negative gcgWeight: pos=%d, weight=%d\n",chromPos,gcgWeight);
			gcgWeight = 100;
		}
		else if (gcgWeight > 5000)
		{
			System.err.printf("Got extra-large gcgWeight: pos=%d, weight=%d\n",chromPos,gcgWeight);
			gcgWeight = 100;
		}
		
		if (hcgWeight < 0)
		{
			System.err.printf("Got negative hcgWeight: pos=%d, weight=%d\n",chromPos,hcgWeight);
			hcgWeight = 100;
		}
		else if (hcgWeight > 5000)
		{
			System.err.printf("Got extra-large hcgWeight: pos=%d, weight=%d\n",chromPos,hcgWeight);
			hcgWeight = 100;
		}
		
		// It's actually divided by two since we're looking at both strands.
		this.gchWeight = (double)gchWeight/2.0;
		this.gcgWeight = (double)gcgWeight/2.0;
		this.hcgWeight = (double)hcgWeight/2.0;
	}
	
	public Cytosine(int chromPos, boolean negStrand, short totalReads, short cReads,
			short cReadsNonconversionFilt, short tReads, short agReads,
			short totalReadsOpposite, short aReadsOpposite, int alleleChromPos,
			String A_BaseUpperCase, String B_BaseUpperCase, short A_CReads, short B_CReads, short A_TReads, short B_TReads, 
			String preBaseRefUpperCase, String nextBaseRefUpperCase ) {
		super();
		this.chromPos = chromPos;
		this.negStrand = negStrand;
		this.totalReads = totalReads;
		this.cReads = cReads;
		this.cReadsNonconversionFilt = cReadsNonconversionFilt;
		this.tReads = tReads;
		this.agReads = agReads;
		this.totalReadsOpposite = totalReadsOpposite;
		this.aReadsOpposite = aReadsOpposite;
		this.alleleChromPos = alleleChromPos;
		this.A_CReads = A_CReads; 
		this.B_CReads = B_CReads;
		this.A_TReads = A_TReads;
		this.B_TReads = B_TReads;
		
		setA_BaseUpperCase(A_BaseUpperCase);
		setB_BaseUpperCase(B_BaseUpperCase);
		setPreBaseRef(preBaseRefUpperCase);
		setNextBaseRef(nextBaseRefUpperCase);

	}
	
	/*** Overridden Comparable methods
	 */

	/* (non-Javadoc)
	 * @see java.lang.Comparable#compareTo(java.lang.Object)
	 */
	public int compareTo(Object o) {

		Cytosine other = (Cytosine)o;
		
		return (this.chromPos - other.chromPos);
	}


	/* (non-Javadoc)
	 * @see java.lang.Object#equals(java.lang.Object)
	 */
	@Override
	public boolean equals(Object obj) {
		return 	(this.compareTo(obj) == 0);
	}
	
	public static Cytosine compositeCytocine(List<Cytosine> inCytocines) throws CloneNotSupportedException
	{
		Cytosine out = null;
		
		for (Cytosine inCytocine : inCytocines)
		{
			if (out == null)
			{
				// First one
				out = (Cytosine)inCytocine.clone();
			}
			else
			{
				out.totalReads += inCytocine.totalReads;
				out.cReads += inCytocine.cReads;
				out.cReadsNonconversionFilt += inCytocine.cReadsNonconversionFilt;
				out.tReads += inCytocine.tReads;
				out.agReads += inCytocine.agReads;
				out.totalReadsOpposite += inCytocine.totalReadsOpposite;
				out.aReadsOpposite += inCytocine.aReadsOpposite;
			}
		
		}
		
		return out;
	}
	
	public double fracNextBaseG()
	{
		if(((double)this.nextBaseTotalReads)!=0){
			return (double)this.nextBaseGreads / ((double)this.nextBaseTotalReads);
		}
		else{
			return -0.01;
		}
		
	}
	
	public double fracPreBaseG()
	{
		if(((double)this.preBaseTotalReads)!=0){
			return (double)this.preBaseGreads / ((double)this.preBaseTotalReads);
		}
		else{
			return -0.01;
		}
		
	}

	public double fracOppositeA()
	{
		return (double)this.aReadsOpposite / ((double)this.aReadsOpposite + (double)this.totalReadsOpposite);
	}
	
	public boolean passesOppositeAFilterDefault()
	{
		return passesOppositeAFilter(DEFAULT_FAIL_MIN_FRAC_A, DEFAULT_FAIL_MIN_ABS_A);
	}
	
	public boolean passesOppositeAFilter(double failMinFracA, int failMinAbsA)
	{
		double fracOppA = this.fracOppositeA();
		boolean failFrac = (!Double.isNaN(fracOppA)) && (fracOppA>=failMinFracA);
		boolean failAbs = (this.aReadsOpposite>=failMinAbsA);
		return !(failFrac && failAbs);
	}

	public double fracMeth(){
		return methyDens;
	}
	
	public double fracMeth(boolean useNonconvFilt)
	{
		int cs = (useNonconvFilt) ? this.cReads : (this.cReads + this.cReadsNonconversionFilt); 
		if(((double)cs + (double)this.tReads) != 0){
			return (double)cs / ((double)cs + (double)this.tReads);
		}
		else{
			return -0.01;
		}
		
	}
	
	public double fracA_Meth()
	{ 
		return (double)this.A_CReads / ((double)this.A_CReads + (double)this.A_TReads);
	}
	
	public double fracB_Meth()
	{ 
		return (double)this.B_CReads / ((double)this.B_CReads + (double)this.B_TReads);
	}
	
	public String variableStepWigLine(boolean useNonconvFilter)
	{
		return String.format("%d\t%f",chromPos, this.fracMeth(useNonconvFilter));
	}

	
	public static void outputCytocinesToFile(PrintWriter pw, Map<Integer, Cytosine> cytocineMap, String prefix, String sampleName, String chr)
	throws IOException
	{
		//System.err.println("About to write " + cpgMap.size() + " Cytosines to file");
		Iterator<Cytosine> cytocineIt = cytocineMap.values().iterator();
		 while (cytocineIt.hasNext())
		{
			Cytosine cytocine = cytocineIt.next();
			//if ((cytocine.getNextBaseRef() != 'G'))// && (cpg.getNextBaseRef() != '0'))
			//{
				
				//	continue CPG;
				
			//}
			//if ((cytocine.getLastBaseRef() == 'G') && (cytocine.getNextBaseRef() != 'G')){
				
			//}
			//else if((cytocine.getLastBaseRef() == 'G') && (cytocine.getNextBaseRef() == 'G')){
				
			//}
			//else if((cytocine.getLastBaseRef() != 'G') && (cytocine.getNextBaseRef() == 'G')){
				
			//}
			//else{
				
			//}
			
			String line = cytocine.toString();
			pw.println(line);
		}
	}
	
	public static void outputCytocinesToFile(PrintWriter pw, Map<String, Cytosine> cytocineMap, String prefix, String sampleName, String chr, boolean asmFlag)
	throws IOException
	{
		//System.err.println("About to write " + cpgMap.size() + " Cytosines to file");
		Iterator<Cytosine> cytocineIt = cytocineMap.values().iterator();
		 while (cytocineIt.hasNext())
		{
			Cytosine cytocine = cytocineIt.next();
			
			String line = cytocine.toString(asmFlag);
			pw.println(line);
		}
	}
	
	public static void outputCytocinesToDb(Map<String, Cytosine> cytocineMap, String tableName, boolean asmFlag)
	throws IOException
	{
		//System.err.println("About to write " + cpgMap.size() + " Cytosines to file");
		Iterator<Cytosine> cytocineIt = cytocineMap.values().iterator();
		 while (cytocineIt.hasNext())
		{
			Cytosine cytocine = cytocineIt.next();
			
			String line = cytocine.toString(asmFlag,true);
			String insertString = "INSERT INTO " + tableName +
								" VALUES (" + line + ")";
			
			try {
				//step 4: create a statement
				Statement stmt = conn.createStatement();
			//	testResult = stmt.executeQuery(queryString);
				
				//step 5: execute a query or update.
				stmt.executeUpdate(insertString);

			}catch(SQLException ex) {
							System.err.println("Insertion: " + ex.getMessage());
			}
			
		}
	}
	
	public static void outputCytocinesToDb(Map<Integer, Cytosine> cytocineMap, String tableName, int minCphCoverage, double minCphMethFrac)
	throws IOException
	{
		//System.err.println("About to write " + cpgMap.size() + " Cytosines to file");
		Iterator<Cytosine> cytocineIt = cytocineMap.values().iterator();
		CYT: while (cytocineIt.hasNext())
		{
			Cytosine cytocine = cytocineIt.next();
			if ((cytocine.getNextBaseRef() != 'G') && (cytocine.getPreBaseRef() != 'G'))// && (cpg.getNextBaseRef() != '0'))
			{
				if ((cytocine.totalReads < minCphCoverage) || (cytocine.fracMeth(true)<minCphMethFrac) || (!cytocine.passesOppositeAFilterDefault()))
				{
					continue CYT;
				}
			}
			String line = cytocine.toStringExpanded(true);
			String insertString = "INSERT INTO " + tableName +
								" VALUES (" + line + ")";
			
			try {
				//step 4: create a statement
				Statement stmt = conn.createStatement();
			//	testResult = stmt.executeQuery(queryString);
				
				//step 5: execute a query or update.
				stmt.executeUpdate(insertString);

			}catch(SQLException ex) {
							System.err.println("Insertion: " + ex.getMessage());
			}
			
		}
	}
	
	public static void outputCytocinesToDb(SortedMap<Integer,TreeMap<Integer,Cytosine>> cytocineMap, String tableName, boolean asmFlag)
	throws IOException
	{
		//System.err.println("About to write " + cpgMap.size() + " Cytosines to file");
		Iterator<TreeMap<Integer,Cytosine>> alleleCytocineIt = cytocineMap.values().iterator();
		 while (alleleCytocineIt.hasNext())
		{
			 TreeMap<Integer,Cytosine> alleleCytocine = alleleCytocineIt.next();
			 Iterator<Cytosine> cytocineIt = alleleCytocine.values().iterator();
			 while(cytocineIt.hasNext()){
				 Cytosine cytocine = cytocineIt.next();
				 String line = cytocine.toString(asmFlag,true);
					String insertString = "INSERT INTO " + tableName +
										" VALUES (" + line + ")";
					
					
					try {
						//step 4: create a statement
						Statement stmt = conn.createStatement();
					//	testResult = stmt.executeQuery(queryString);
						
						//step 5: execute a query or update.
						stmt.executeUpdate(insertString);

					}catch(SQLException ex) {
									System.err.println("Insertion: " + ex.getMessage());
					}
				 
			 }
			
			
		}
	}
	
	public static PrintWriter outputChromToFile(Map<Integer,Cytosine> cytocineMap, String prefix, String sampleName, String chr)
	throws IOException
	{
		
		String fn = prefix + sampleName + "_" + chr + ".txt";
		PrintWriter writer = new PrintWriter(new File(fn));
		
		outputCytocinesToFile(writer, cytocineMap, prefix, sampleName, chr);
		
		//writer.close();
		return writer;
	}
	
	
	public static PrintWriter outputChromToFile(Map<String, Cytosine> cytocineMap, String prefix, String sampleName, String chr, boolean asmFlag)
	throws IOException
	{
		
		String fn = prefix + sampleName + "_" + chr + ".txt";
		PrintWriter writer = new PrintWriter(new File(fn));
		
		outputCytocinesToFile(writer, cytocineMap, prefix, sampleName, chr, asmFlag);
		
		//writer.close();
		return writer;
	}

	public static void outputChromToDb(Map<String, Cytosine> cytocineMap, String tableName, Connection cConn, boolean asmFlag)
	throws IOException
	{
		
		conn=cConn;
		creatTableInDb(tableName, asmFlag);
		outputCytocinesToDb(cytocineMap, tableName, asmFlag);

	}
	
	public static void outputChromToDb(Map<Integer, Cytosine> cytocineMap, String tableName, Connection cConn, int minCphCoverage, double minCphMethFrac)
	throws IOException
	{
		
		conn=cConn;
		creatTableInDb(tableName);
		outputCytocinesToDb(cytocineMap, tableName, minCphCoverage, minCphMethFrac);

	}
	
	
	
	public static void outputChromToDb(SortedMap<Integer,TreeMap<Integer,Cytosine>> cytocineMap, String tableName, Connection cConn, boolean asmFlag)
	throws IOException
	{
		
		conn=cConn;
		creatTableInDb(tableName, asmFlag);
		outputCytocinesToDb(cytocineMap, tableName, asmFlag);

	}
	
	public static void creatTableInDb(String tableName, boolean asmFlag){
		Statement stmtDrop, stmtCreate, stmtCheck;
		//ResultSet tempCheck;
		String dropString = "DROP TABLE " + tableName;
		String checkString = "SHOW TABLES LIKE " + "\"" + tableName + "\"";
		String createString = "CREATE TABLE " + tableName +
        					"(`chromPos` INT UNSIGNED NOT NULL," +
        					"`strand` enum('+','-') NOT NULL," +
        					"`totalReads` SMALLINT UNSIGNED NOT NULL," +
        					"`cReads` SMALLINT UNSIGNED NOT NULL," +
        					"`cReadsNonconversionFilt` SMALLINT UNSIGNED NOT NULL," +
        					"`tReads` SMALLINT UNSIGNED NOT NULL," +
        					"`agReads` SMALLINT UNSIGNED NOT NULL," +
        					"`totalReadsOpposite` SMALLINT UNSIGNED NOT NULL," +
        					"`aReadsOpposite` SMALLINT UNSIGNED NOT NULL," +
        					"`alleleChromPos` INT NOT NULL," +
        					"`ABaseRefUpperCase` CHAR(1)," +
        					"`BBaseRefUpperCase` CHAR(1)," +
        					"`ACReads` SMALLINT UNSIGNED NOT NULL," +
        					"`BCReads` SMALLINT UNSIGNED NOT NULL," +
        					"`ATReads` SMALLINT UNSIGNED NOT NULL," +
        					"`BTReads` SMALLINT UNSIGNED NOT NULL," +
        					"`preBaseGreads` SMALLINT UNSIGNED NOT NULL default '0'," +
        					"`preBaseTotalReads` SMALLINT UNSIGNED NOT NULL default '0'," +
        					"`preBaseRefUpperCase` CHAR(1)," +
        					"`nextBaseGreads` SMALLINT UNSIGNED NULL default '0'," +
        					"`nextBaseTotalReads` SMALLINT UNSIGNED NULL default '0'," +
        					"`nextBaseRefUpperCase` CHAR(1)," +
        					/*"`fracMeth` FLOAT(5,2) NOT NULL," +
        					"`fracAMeth` FLOAT(5,2) NOT NULL," +
        					"`fracBMeth` FLOAT(5,2) NOT NULL," +
        					"`fracPreBaseG` FLOAT(5,2) NOT NULL," +
        					"`fracNextBaseG` FLOAT(5,2) NOT NULL," +
        					"`gchWeight` FLOAT(5,2) NOT NULL," +
        					"`gcgWeight` FLOAT(5,2) NOT NULL," +
        					"`hcgWeight` FLOAT(5,2) NOT NULL," +*/
        					"PRIMARY KEY (chromPos,alleleChromPos))";
		//String dropString = "DROP TABLE "+ tableName;
		try {
//step 4: create a statement
			stmtCheck = conn.createStatement();
//step 5: execute a query or update.
			//stmt.execute(dropString);//if exists, drop it, get new one
			//tempCheck = stmtCheck.executeQuery(checkString);
			if (stmtCheck.executeQuery(checkString).next()){
				stmtDrop = conn.createStatement();
				stmtDrop.execute(dropString);
				System.err.println( tableName + " exist! Drop the table!");
			}
			stmtCreate = conn.createStatement();
			stmtCreate.executeUpdate(createString);
			System.err.println("CreateTable: " + tableName);

		}catch(SQLException ex) {
			System.err.println("CreateTable: " + ex.getMessage());
		}
	}
	
	public static void creatTableInDb(String tableName){
		Statement stmtDrop, stmtCreate, stmtCheck;
		//ResultSet tempCheck;
		String dropString = "DROP TABLE " + tableName;
		String checkString = "SHOW TABLES LIKE " + "\"" + tableName + "\"";
		String createString = "CREATE TABLE " + tableName +
        					"(`chromPos` INT UNSIGNED NOT NULL," +
        					"`strand` enum('+','-') NOT NULL," +
        					"`totalReads` SMALLINT UNSIGNED NOT NULL," +
        					"`cReads` SMALLINT UNSIGNED NOT NULL," +
        					"`cReadsNonconversionFilt` SMALLINT UNSIGNED NOT NULL," +
        					"`tReads` SMALLINT UNSIGNED NOT NULL," +
        					"`agReads` SMALLINT UNSIGNED NOT NULL," +
        					"`totalReadsOpposite` SMALLINT UNSIGNED NOT NULL," +
        					"`aReadsOpposite` SMALLINT UNSIGNED NOT NULL," +
        					"`preBaseGreads` SMALLINT UNSIGNED NOT NULL default '0'," +
        					"`preBaseTotalReads` SMALLINT UNSIGNED NOT NULL default '0'," +
        					"`preBaseRefUpperCase` CHAR(1)," +
        					"`nextBaseGreads` SMALLINT UNSIGNED NULL default '0'," +
        					"`nextBaseTotalReads` SMALLINT UNSIGNED NULL default '0'," +
        					"`nextBaseRefUpperCase` CHAR(1)," +
        					"`fracMeth` FLOAT(5,2) NOT NULL," +
        					"`fracPreBaseG` FLOAT(5,2) NOT NULL," +
        					"`fracNextBaseG` FLOAT(5,2) NOT NULL," +
        					"`gchWeight` FLOAT(5,2) NOT NULL," +
        					"`gcgWeight` FLOAT(5,2) NOT NULL," +
        					"`hcgWeight` FLOAT(5,2) NOT NULL," +
        					"PRIMARY KEY (chromPos))";
		//String dropString = "DROP TABLE "+ tableName;
		try {
//step 4: create a statement
			stmtCheck = conn.createStatement();
//step 5: execute a query or update.
			//stmt.execute(dropString);//if exists, drop it, get new one
			//tempCheck = stmtCheck.executeQuery(checkString);
			if (stmtCheck.executeQuery(checkString).next()){
				stmtDrop = conn.createStatement();
				stmtDrop.execute(dropString);
				System.err.println( tableName + " exist! Drop the table!");
			}
			stmtCreate = conn.createStatement();
			stmtCreate.executeUpdate(createString);
			System.err.println("CreateTable: " + tableName);

		}catch(SQLException ex) {
			System.err.println("CreateTable: " + ex.getMessage());
		}
	}
	
	@Override
	public String toString() {

		if (true)
		{
			return this.toStringExpanded();
		}
		else
		{
		return String.format("%d\t%c\t%d\t%d\t%d\t%d\t%d\t%d\t%d", 
				chromPos,
				(negStrand) ? '-' : '+',
				totalReads,
				cReads,
				cReadsNonconversionFilt,
				tReads,
				agReads,
				totalReadsOpposite,
				aReadsOpposite
				);
		}
	}
	
	public String toString(boolean asmFlag) {

		return String.format("%d\t%c\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%c\t%c\t%d\t%d\t%d\t%d\t%d\t%d\t%c\t%d\t%d\t%c\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f", 
				chromPos,
				(negStrand) ? '-' : '+',
				totalReads,
				cReads,
				cReadsNonconversionFilt,
				tReads,
				agReads,
				totalReadsOpposite,
				aReadsOpposite,
				alleleChromPos, 
				A_BaseUpperCase, 
				B_BaseUpperCase, 
				A_CReads,
				B_CReads,
				A_TReads,
				B_TReads,
				preBaseGreads,
				preBaseTotalReads,
				preBaseRefUpperCase,
				nextBaseGreads,
				nextBaseTotalReads,
				nextBaseRefUpperCase,
				100*this.fracMeth(true),
				100*this.fracA_Meth(),
				100*this.fracB_Meth(),
				100*this.fracPreBaseG(),
				100*this.fracNextBaseG(),
				this.gchWeight,
				this.gcgWeight,
				this.hcgWeight
				);
	}
	
	public String toString(boolean asmFlag, boolean online) {

		return String.format("%d,'%c',%d,%d,%d,%d,%d,%d,%d,%d,'%c','%c',%d,%d,%d,%d,%d,%d,'%c',%d,%d,'%c'", 
				chromPos,
				(negStrand) ? '-' : '+',
				totalReads,
				cReads,
				cReadsNonconversionFilt,
				tReads,
				agReads,
				totalReadsOpposite,
				aReadsOpposite,
				alleleChromPos, 
				A_BaseUpperCase, 
				B_BaseUpperCase, 
				A_CReads,
				B_CReads,
				A_TReads,
				B_TReads,
				preBaseGreads,
				preBaseTotalReads,
				preBaseRefUpperCase,
				nextBaseGreads,
				nextBaseTotalReads,
				nextBaseRefUpperCase
			/*	100*this.fracMeth(true),
				100*this.fracA_Meth(),
				100*this.fracB_Meth(),
				100*this.fracPreBaseG(),
				100*this.fracNextBaseG(),
				this.gchWeight,
				this.gcgWeight,
				this.hcgWeight*/
				);
	}
	
	public String toStringExpanded(boolean online) 
	{
			return String.format("%d,'%c',%d,%d,%d,%d,%d,%d,%d,%d,%d,'%c',%d,%d,'%c',%.2f,%.2f,%.2f,%.2f,%.2f,%.2f", 
					chromPos,
					(negStrand) ? '-' : '+',
					totalReads,
					cReads,
					cReadsNonconversionFilt,
					tReads,
					agReads,
					totalReadsOpposite,
					aReadsOpposite,
					preBaseGreads,
					preBaseTotalReads,
					preBaseRefUpperCase,
					nextBaseGreads,
					nextBaseTotalReads,
					nextBaseRefUpperCase,
					100*this.fracMeth(true),
					100*this.fracPreBaseG(),
					100*this.fracNextBaseG(),
					this.gchWeight,
					this.gcgWeight,
					this.hcgWeight
					);
			
	}
	
	public String toStringExpanded() 
	{

		return String.format("%d\t%c\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%c\t%d\t%d\t%c\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f", 
				chromPos,
				(negStrand) ? '-' : '+',
				totalReads,
				cReads,
				cReadsNonconversionFilt,
				tReads,
				agReads,
				totalReadsOpposite,
				aReadsOpposite,
				preBaseGreads,
				preBaseTotalReads,
				preBaseRefUpperCase,
				nextBaseGreads,
				nextBaseTotalReads,
				nextBaseRefUpperCase,
				100*this.fracMeth(true),
				100*this.fracPreBaseG(),
				100*this.fracNextBaseG(),
				this.gchWeight,
				this.gcgWeight,
				this.hcgWeight
				);
				
	}
	

	@Override
	protected Object clone() throws CloneNotSupportedException {
		Cytosine out = (Cytosine)super.clone();
		
		System.err.printf("Cloning cytocine, old tReads=%d, new tReads=%d\n",  this.tReads, out.tReads);
		return out;
	}

	public double getGchWeight() {
		return gchWeight;
	}
	
	public double getGcgWeight() {
		return gcgWeight;
	}
	
	public double getHcgWeight() {
		return hcgWeight;
	}

//	public void setCpgWeight(double cpgWeight) {
//		this.cpgWeight = cpgWeight;
//	}

	public char getNextBaseRef() {
		return nextBaseRefUpperCase;
	}

	public void setNextBaseRef(char nextBaseRef) {
		this.nextBaseRefUpperCase = Character.toUpperCase(nextBaseRef);
	}
	
	public void setPreBaseRef(String preBaseRef){
		this.preBaseRefUpperCase = preBaseRef.charAt(0);
	}
	
	public void setNextBaseRef(String nextBaseRef){
		this.nextBaseRefUpperCase = nextBaseRef.charAt(0);
	}
	
	public char getPreBaseRef() {
		return preBaseRefUpperCase;
	}

	public void setPreBaseRef(char preBaseRef) {
		this.preBaseRefUpperCase = Character.toUpperCase(preBaseRef);
	}

	public char getA_BaseUpperCase() {
		return A_BaseUpperCase;
	}
	
	public void setA_BaseUpperCase(String A_BaseUpperCase){
		this.A_BaseUpperCase = A_BaseUpperCase.charAt(0);
	}
	
	public void setA_BaseUpperCase(char A_BaseUpperCase){
		this.A_BaseUpperCase = Character.toUpperCase(A_BaseUpperCase);
	}
	
	public char getB_BaseUpperCase() {
		return B_BaseUpperCase;
	}
	
	public void setB_BaseUpperCase(String B_BaseUpperCase){
		this.B_BaseUpperCase = B_BaseUpperCase.charAt(0);
	}
	
	public void setB_BaseUpperCase(char B_BaseUpperCase){
		this.B_BaseUpperCase = Character.toUpperCase(B_BaseUpperCase);
	}
	
	public short getA_CReads() {
		return this.A_CReads;
	}
	
	public short getA_TReads() {
		return this.A_TReads;
	}
	
	public short getB_CReads() {
		return this.B_CReads;
	}
	
	public short getB_TReads() {
		return this.B_TReads;
	}
	
	public StrandedFeature.Strand getStrand()
	{
		return (this.negStrand) ? StrandedFeature.NEGATIVE : StrandedFeature.POSITIVE;
	}
	
	public String getStrandStr()
	{
		return (this.negStrand) ? "-" : "+";
	}
	
	public short totalReadsC(boolean useNonconvFilt)
	{
		int cs = (useNonconvFilt) ? this.cReads : (this.cReads + this.cReadsNonconversionFilt); 
		return (short)cs;
	}
	
	public short totalReadsCorT(boolean useNonconvFilt)
	{
		short cs = totalReadsC(useNonconvFilt);
		return (short)(cs + this.tReads);
	}

	public Cytosine downsample(double downsamplingFactor) {
		// TODO Auto-generated method stub
		
		System.err.println("Cytocine::downsample not yet supported");
		System.exit(1);
		
		return null;
	}
	
	
	
	
}
