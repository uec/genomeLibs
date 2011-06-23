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
import org.biojava.bio.seq.StrandedFeature;
import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;
import org.usckeck.genome.ChromFeatures;
import org.usckeck.genome.GFFUtils;

import com.sun.org.apache.xalan.internal.xsltc.compiler.Pattern;

import edu.usc.epigenome.genomeLibs.BiojavaUtils;
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
	@Option(name="-cpgMethTableFormat",usage="If used, this creates tables of the MethylCGRich format , one for each feature type and chrom (default false)")
    protected boolean cpgMethTableFormat = false;
	@Option(name="-methTableUseMidpoint",usage="Makes a single entry at the midpoint of each element (default is to put one for start point and one for end)")
    protected boolean methTableUseMidpoint = false;

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
				// print the list of available options
				parser.printUsage(System.err);
				System.err.println();
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
			
			// If we're doing meth tables, it's one per feat type, so we can close output files
			if (this.cpgMethTableFormat) this.closeOutputFiles();
		}
	
		// Finish up
		this.closeOutputFiles();
	}
		
	protected void closeOutputFiles()
	{
		// Close files
		for (PrintWriter writer : this.outFiles.values())
		{
			writer.close();
		}
	}
	
	public void outputFeat(GFFRecord feat, String featType)
	throws IOException
	{
		if (this.cpgMethTableFormat)
		{
			outputFeatMethTable(feat,featType);
		}
		else
		{
			outputFeatFeatTable(feat,featType);
		}
			
	}
	
	public void outputFeatMethTable(GFFRecord feat, String featType)
	throws IOException
	{
		String chr = feat.getSeqName();
		
		// If it's oriented, we make a separate table for starts and ends.  Unless
		// it's trivial size (<2bp), in which case we only put it into starts
		int featLeft = feat.getStart();
		int featRight = feat.getEnd();
		int featLen = featRight-featLeft+1;
		boolean featRev = (feat.getStrand()==StrandedFeature.NEGATIVE);

		
		if (this.methTableUseMidpoint)
		{
			int midpoint = Math.round(((float)featLeft+(float)featRight)/2.0f);
			outputFeatMethTableEnd(chr, featType, midpoint, featLen, "Midpoint", feat.getStrand());
		}
		else if (featLen <= 10)
		{
			// A single point basically
			outputFeatMethTableEnd(chr, featType, (featRev) ? featRight : featLeft, featLen, "Starts", feat.getStrand());
		}
		else if (feat.getStrand() == StrandedFeature.UNKNOWN)
		{
			// Unstranded
			outputFeatMethTableEnd(chr, featType, featLeft, featLen, "Unstranded", StrandedFeature.POSITIVE);
			outputFeatMethTableEnd(chr, featType, featRight, featLen, "Unstranded", StrandedFeature.NEGATIVE);
		}
		else
		{
			// Stranded
			outputFeatMethTableEnd(chr, featType, featLeft, featLen, (featRev)?"Ends":"Starts", StrandedFeature.POSITIVE);
			outputFeatMethTableEnd(chr, featType, featRight, featLen, (featRev)?"Starts":"Ends", StrandedFeature.NEGATIVE);
		}
	
	}
	
	public void outputFeatMethTableEnd(String chr, String featType, int chromPos, int featLen, String endLabel, StrandedFeature.Strand strand)
	throws IOException
	{
		String fn = String.format("methylCGsRich_%s_%s_%s.txt", featType, endLabel, chr);
		
		PrintWriter writer = this.outFiles.get(fn);
		if (writer == null)
		{
			writer = new PrintWriter(new File(fn));
			this.outFiles.put(fn, writer);
		}
		
		
		List<String> flds = new ArrayList<String>(12);
		
		flds.add(Integer.toString(chromPos));
		flds.add(Character.toString(BiojavaUtils.strandToSymbol(strand)));
		flds.add("1");
		flds.add("1");
		flds.add("0");
		flds.add("0");
		flds.add("0");
		flds.add("0");
		flds.add("0");
		flds.add("0");
		flds.add("0");
		flds.add(Integer.toString(featLen));
		
		writer.println( ListUtils.excelLine(flds.toArray(new String[1])));
	}
	
	public void outputFeatFeatTable(GFFRecord feat, String featType)
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
		
		String refseqId = "NULL";
		String refseqVers = "NULL";
		String name = "NULL";
		String featName = GFFUtils.getGffRecordName(feat);

		if (featName != null) name = featName;
	
		flds.add(refseqId);
		flds.add(refseqVers);
		flds.add(name);
		
		writer.println( ListUtils.excelLine(flds.toArray(new String[1])));
	}
	
//	public static void outputChromToFile(Map<Integer,Cpg> cpgMap, String prefix, String sampleName, String chr)
//	throws IOException
//	{
//		
//		String fn = prefix + sampleName + "_" + chr + ".txt";
//		PrintWriter writer = new PrintWriter(new File(fn));
//		
//		Iterator<Cpg> cpgIt = cpgMap.values().iterator();
//		while (cpgIt.hasNext())
//		{
//			Cpg cpg = cpgIt.next();
//			String line = cpg.toString();
//			writer.println(line);
//		}
//		writer.close();
//	}
}
