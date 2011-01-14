package edu.usc.epigenome.scripts;

import java.io.File;
import java.io.FileOutputStream;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;

import org.biojava.bio.program.gff.GFFRecord;
import org.biojava.bio.program.gff.SimpleGFFRecord;
import org.biojava.bio.seq.StrandedFeature;

import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;
import org.usckeck.genome.ChromFeatures;

import sun.tools.tree.ThisExpression;

import edu.usc.epigenome.genomeLibs.GFFUtils;
import edu.usc.epigenome.genomeLibs.GoldAssembly;
import edu.usc.epigenome.genomeLibs.ListUtils;
import edu.usc.epigenome.genomeLibs.MatUtils;
import edu.usc.epigenome.genomeLibs.ChromScores.ChromScoresArrayInt;
import edu.usc.epigenome.genomeLibs.FeatAligners.FeatAligner;
import edu.usc.epigenome.genomeLibs.FeatAligners.FeatAlignerAveraging;
import edu.usc.epigenome.genomeLibs.FeatAligners.FeatAlignerEachfeat;
import edu.usc.epigenome.genomeLibs.FeatDb.FeatDbQuerier;
import edu.usc.epigenome.genomeLibs.FeatDb.FeatIterator;
import edu.usc.epigenome.genomeLibs.MethylDb.Cpg;
import edu.usc.epigenome.genomeLibs.MethylDb.CpgIterator;
import edu.usc.epigenome.genomeLibs.MethylDb.CpgIteratorMultisample;
import edu.usc.epigenome.genomeLibs.MethylDb.MethylDbQuerier;
import edu.usc.epigenome.genomeLibs.MethylDb.MethylDbUtils;


public class FeatCounter {

	private static final String C_USAGE = "Use: FeatCounter elements.gtf";
	
//	@Option(name="-skipUnoriented",usage="If set, skip any unoriented feature (default false)")
//	protected boolean skipUnoriented = false;
	@Option(name="-chrom",multiValued=true,usage="One or more chroms, eg. -chrom chr1 -chrom chr5")
	protected List<String> chrs = new ArrayList<String>(25);
	@Option(name="-feature",multiValued=true,usage="One or more features from features_ tables")
	protected List<String> features = new ArrayList<String>(25);
	@Option(name="-flank",multiValued=false,usage="Flanking sequence to search for each feature (default 0=exact overlap)")
	protected int flank = 0;
	@Option(name="-includeSummaryCounts",multiValued=false,usage="If set, the last line are total counts for each column")
	protected boolean includeSummaryCounts = false;
	@Option(name="-tssMaxFlank",multiValued=false,usage="This is the maximum flank we will use to get TSS expression")
	protected int tssMaxFlank = 500000;
	@Option(name="-expressionTerm",usage="One or more expression from the infiniumExpr_chr table")
	protected List<String> expressionTerms = new ArrayList<String>(5);
	@Option(name="-includeTssAndBorrowExpression",usage="If set, tss is the first feat, and we get expression from there")
	protected boolean includeTssAndBorrowExpression = false;

	// receives other command line parameters than options
	@Argument
	private List<String> arguments = new ArrayList<String>();


	// class vars


	/**
	 * @param args
	 */
	public static void main(String[] args)
	throws Exception	
	{
		Logger.getLogger(Logger.GLOBAL_LOGGER_NAME).setLevel(Level.INFO);
		new FeatCounter().doMain(args);
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

		if (arguments.size()>1)
		{
			System.err.println(C_USAGE);
			parser.printUsage(System.err);
			return;
		}
		
		if (features.size()<1)
		{
			System.err.println("Must supply at least 1 -feature argument");
			System.err.println(C_USAGE);
			parser.printUsage(System.err);
			return;
		}		
		
		if (chrs.size()==0) chrs = MethylDbUtils.CHROMS;
		
		int tssFlank = Math.max(this.tssMaxFlank, this.flank);
		
		ChromFeatures targetFeats = null;
		if (arguments.size()>=1)
		{
			String targetFn = arguments.get(0);
			targetFeats = new ChromFeatures(targetFn, true);
		}

		// Go through target feats one by one
		int[] totals = new int[features.size()+1];
		ListUtils.setDelim(",");
		String expressionSec = (this.expressionTerms.size()==0) ? "" : ("," + ListUtils.excelLine(this.expressionTerms));
		String tssSec = (this.includeTssAndBorrowExpression) ? ",tss" : "";
		System.out.printf("%s,%s,%s%s%s,%s,%s\n","chrom","start","end",expressionSec, tssSec, ListUtils.excelLine(features),"No overlap");
		for (String chrStr : chrs)
		{
			// Some of my files have these forms.
			if (chrStr.equalsIgnoreCase("chrX")) chrStr = "chrX";
			if (chrStr.equalsIgnoreCase("chry")) chrStr = "chrY";
			if (chrStr.equalsIgnoreCase("chrm")) chrStr = "chrM";

			if (targetFeats == null)
			{
				int chrLen = GoldAssembly.chromLengthStatic(chrStr, "hg18");
				
				// Take size of features.
				for (int n = 0; n < features.size(); n++)
				{
					ChromScoresArrayInt scores = new ChromScoresArrayInt("hg18");
					String featType = features.get(n);

					// Get the overlapping feats
					FeatDbQuerier params = new FeatDbQuerier();
					params.addRangeFilter(chrStr, 0, 0); // Set to 0,0 for whole chromosome
					params.addFeatFilter(featType);
					FeatIterator feats = new FeatIterator(params);
					while (feats.hasNext())
					{
						GFFRecord rec = feats.next();
						scores.addRange(chrStr, Math.max(0,rec.getStart()-flank), Math.min(chrLen-1,rec.getEnd()+flank), 1);
					}

					totals[n] += scores.getCoverageTotal(chrStr);
				}

			}
			else
			{			
				int nF = 0;
				int chrNum = (new ChromFeatures()).chrom_from_public_str(chrStr);
				Iterator targetit = targetFeats.featureIterator(chrNum);
				while (targetit.hasNext())
				{

					SimpleGFFRecord target = (SimpleGFFRecord) targetit.next();
					//System.err.printf("Working on target: %s\n",GFFUtils.gffBetterString(target));

					boolean overlapsSomething = false;
					int start = (includeTssAndBorrowExpression) ? -1 : 0; 
					for (int n = start; n < features.size(); n++)
					{
						String featType = (n>=0) ? features.get(n) : "tss";

						// Get the overlapping feats
						FeatDbQuerier params = new FeatDbQuerier();
						int thisFlank = (n>=0) ? flank : tssFlank;
						params.addRangeFilter(target.getSeqName(), target.getStart()-thisFlank, target.getEnd()+thisFlank);
						params.addFeatFilter(featType);
						FeatIterator feats = new FeatIterator(params);

						// Find the closest feature
						GFFRecord closestFeat = null;
						int closestFeatDist = Integer.MAX_VALUE;
						int featsSearched = 0;
						while (feats.hasNext())
						{
							GFFRecord feat = feats.next(); // This is the "subject" feat (like tss)
							StrandedFeature.Strand featStrand = feat.getStrand();
							if (featStrand == StrandedFeature.UNKNOWN)
							{
//								System.err.printf("Why does feat have no strand? %s\n", GFFUtils.gffBetterString(feat));
							}
							
							int dist;
							if (feat.getStart() > target.getEnd())
							{
								dist = feat.getStart() - target.getEnd();
								if (featStrand == StrandedFeature.POSITIVE) dist *= -1;
							}
							else if (feat.getEnd() < target.getStart())
							{
								dist = target.getStart() - feat.getEnd();
								if (featStrand == StrandedFeature.NEGATIVE) dist *= -1;
							}
							else
							{
								// It's within the bounds
								dist = 0;
							}
							if (Math.abs(dist)<Math.abs(closestFeatDist))
							{
								closestFeatDist = dist;
								closestFeat = feat;
							}
							featsSearched++;
//							System.err.printf("\tEvaluating feat with dist=%d (closest=%s): %s\n",dist, (closestFeat==feat),GFFUtils.gffBetterString(feat));
						}
//						System.err.printf("Searched %d overlapping %s feats for: %s\n", featsSearched, featType, GFFUtils.gffBetterString(target));
						
						
//						boolean overlap = (feats.hasNext());
						boolean overlap = (Math.abs(closestFeatDist) <= this.flank);
						
						
						if ((n>=0) && overlap) totals[n]++;
						overlapsSomething |= overlap;

						// Get the expression.  If we are using TSS, we have to make sure we have one within distance
						if (((n==0) && !this.includeTssAndBorrowExpression) || ((n==-1) && this.includeTssAndBorrowExpression))
						{
							System.out.printf("%d,%d,%d",chrNum,target.getStart(),target.getEnd());

							GFFRecord exprTarget = (this.includeTssAndBorrowExpression) ? closestFeat : target;
							
							for (String expressionTerm : expressionTerms)
							{
								double term = Double.NaN; 
								if (exprTarget != null) term = MethylDbUtils.fetchMeanExpression(chrStr, GFFUtils.getGffRecordName(exprTarget),expressionTerm);
								//System.out.printf(",%.5f (n=%d, name=%s)", term, n, GFFUtils.getGffRecordName(exprTarget));
								System.out.printf(",%.5f", term);
							}
						}
						
						// Now print the overlap. 
						//System.out.printf(",%d",(overlap)?1:0);  // The old way, boolean
						if (closestFeat == null)
						{
							System.out.printf(",NaN");
						}
						else
						{
							System.out.printf(",%d",closestFeatDist);
						}
						
						
					}
					System.out.printf(",%d",(!overlapsSomething)?1:0);
					if (!overlapsSomething) totals[features.size()]++;
					System.out.println();
				}
			}
		}
		
		ListUtils.setDelim(",");
		if (includeSummaryCounts)
		{
			System.out.printf("%s,%s,%s,%s\n","totals","a","b",ListUtils.excelLine(totals));
		}
	}
	
	
}
