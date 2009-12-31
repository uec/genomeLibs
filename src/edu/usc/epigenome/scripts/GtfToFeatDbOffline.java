package edu.usc.epigenome.scripts;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;

import org.biojava.bio.program.gff.GFFRecord;
import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;
import org.usckeck.genome.ChromFeatures;
import org.usckeck.genome.GFFUtils;

import com.sun.org.apache.xalan.internal.xsltc.compiler.Pattern;

import edu.usc.epigenome.genomeLibs.ListUtils;
import edu.usc.epigenome.genomeLibs.FeatDb.FeatDbQuerier;
import edu.usc.epigenome.genomeLibs.MethylDb.Cpg;
import edu.usc.epigenome.genomeLibs.MethylDb.CpgIterator;
import edu.usc.epigenome.genomeLibs.MethylDb.MethylDbQuerier;



public class GtfToFeatDbOffline {

	private static final String C_USAGE = "Use: GtfToFeatDb -tablePrefix " + FeatDbQuerier.DEFAULT_TABLE_PREFIX + 
	" feat1.gtf file2.gtf ...";
	
	@Option(name="-tablePrefix",usage="Prefix for DB table (default " + FeatDbQuerier.DEFAULT_TABLE_PREFIX + ")")
    protected String tablePrefix = FeatDbQuerier.DEFAULT_TABLE_PREFIX;
 	// receives other command line parameters than options
	@Argument
	private List<String> arguments = new ArrayList<String>();

	// Object vars
	Map<String,PrintWriter> outFiles = new HashMap<String,PrintWriter>(30);

	
	public static void main(String[] args)
	throws Exception	
	{
		new GtfToFeatDbOffline().doMain(args);
	}

	public void doMain(String[] args)
	throws Exception
	{
		CmdLineParser parser = new CmdLineParser(this);
		// if you have a wider console, you could increase the value;
		// here 80 is also the default
		int chrSt = 0, chrEnd = 0;
		String chr;
		parser.setUsageWidth(80);
		try
		{
			parser.parseArgument(args);

			if(arguments.size() < 1 ) {
				System.err.println(C_USAGE);
				System.exit(1);
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
		
		Logger.getLogger(Logger.GLOBAL_LOGGER_NAME).setLevel(Level.SEVERE);
		
		// Setup for file writing
		ListUtils.setDelim("\t");

		
		// Go through input files
		int onFeat = 0;
		int nF = arguments.size(); 
		for (String fn : arguments)
		{
			String fnShort = (new File(fn)).getName();
			String featType = fnShort.substring(0,Math.min(39, fnShort.length()));
			Logger.getLogger(Logger.GLOBAL_LOGGER_NAME).severe(String.format("On feat %d/%d: %s\n", ++onFeat, nF,featType));
			
			ChromFeatures feats = new ChromFeatures(fn, true);
			Iterator featit = feats.featureIterator();
			while (featit.hasNext())
			{
				GFFRecord rec = (GFFRecord) featit.next();
				this.outputFeat(rec, featType);
			}
		}
	
		// Close files
		for (PrintWriter writer : this.outFiles.values())
		{
			writer.close();
		}
		
	}
		
	public void outputFeat(GFFRecord feat, String featType)
	throws IOException
	{
		String chr = feat.getSeqName();
		String fn = this.tablePrefix +  chr + ".txt";
		
		PrintWriter writer = this.outFiles.get(fn);
		if (writer == null)
		{
			writer = new PrintWriter(new File(fn));
			this.outFiles.put(fn, writer);
		}
		
		List<String> flds = new ArrayList<String>(10);
		
		flds.add(featType);
		flds.add(Integer.toString(feat.getStart()));
		flds.add(Integer.toString(feat.getEnd()));
		flds.add(Character.toString(feat.getStrand().getToken()));

		double score = feat.getScore();
		flds.add( (Double.isInfinite(score) || Double.isNaN(score)) ? "." : Double.toString(score));
		
		String featName = GFFUtils.getGffRecordName(feat);
		if (featName==null) featName = "NULL";
		flds.add(featName);
		
		writer.println( ListUtils.excelLine(flds.toArray(new String[1])));
	}
	
	public static void outputChromToFile(Map<Integer,Cpg> cpgMap, String prefix, String sampleName, String chr)
	throws IOException
	{
		
		String fn = prefix + sampleName + "_" + chr + ".txt";
		PrintWriter writer = new PrintWriter(new File(fn));
		
		Iterator<Cpg> cpgIt = cpgMap.values().iterator();
		while (cpgIt.hasNext())
		{
			Cpg cpg = cpgIt.next();
			String line = cpg.toString();
			writer.println(line);
		}
		writer.close();
	}
}
