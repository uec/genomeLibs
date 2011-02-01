package edu.usc.epigenome.scripts;

import java.io.File;
import java.io.FileOutputStream;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashSet;
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
import edu.usc.epigenome.genomeLibs.ListUtils;
import edu.usc.epigenome.genomeLibs.MatUtils;
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


public class FeatDbReorientGtfFeats {

	private static final String C_USAGE = "Use: FeatDbReorientGtfFeats -flank 0 -chrom chr1 -chrom chr2 featsToBeReoriented.gtf referenceFeatType";
	
//	@Option(name="-skipUnoriented",usage="If set, skip any unoriented feature (default false)")
//	protected boolean skipUnoriented = false;
	@Option(name="-chrom",multiValued=true,usage="One or more chroms, eg. -chrom chr1 -chrom chr5")
	protected List<String> chrs = new ArrayList<String>(25);
	@Option(name="-flank",multiValued=true,usage="Search this many base pairs flanking the features of interes")
	protected int flank = 0;
	@Option(name="-dontRename",usage="If unset, we rename elements to match the feature used to reorient")
	protected boolean dontRename = false;

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
		new FeatDbReorientGtfFeats().doMain(args);
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

		if (arguments.size()!=2)
		{
			System.err.println(C_USAGE);
			parser.printUsage(System.err);
			return;
		}
		
		if (chrs.size()==0) chrs = MethylDbUtils.CHROMS;
		String targetFn = arguments.get(0);
		String featType = arguments.get(1);
		
		String outFn = (new File(targetFn)).getName().replaceFirst(".[^.]*$", "") + ".reorientedBy." +
		((featType.length()>15) ? featType.substring(0,15) : featType) + ".gtf";
		Logger.getAnonymousLogger().severe(String.format("%s,%s => %s\n", targetFn,featType,outFn));

		ChromFeatures outFeats = new ChromFeatures();
		
		// Go through target feats one by one
		ChromFeatures targetFeats = new ChromFeatures(targetFn, true);
		for (String chrStr : chrs)
		{
			// Some of my files have these forms.
			if (chrStr.equalsIgnoreCase("chrX")) chrStr = "chrX";
			if (chrStr.equalsIgnoreCase("chry")) chrStr = "chrY";
			if (chrStr.equalsIgnoreCase("chrm")) chrStr = "chrM";
			
			Iterator targetit = targetFeats.featureIterator((new ChromFeatures()).chrom_from_public_str(chrStr));
			while (targetit.hasNext())
			{
				SimpleGFFRecord target = (SimpleGFFRecord) targetit.next();

				// Get the overlapping feats
				FeatDbQuerier params = new FeatDbQuerier();
				params.addRangeFilter(target.getSeqName(), target.getStart()-this.flank, target.getEnd()+this.flank);
				params.addFeatFilter(featType);
				FeatIterator feats = new FeatIterator(params);

				// Figure out the strand
				boolean overlapsFw = false;
				boolean overlapsRev = false;
				int nOvs = 0;
				HashSet<String> newNames = new HashSet<String>(20);
				while (feats.hasNext())
				{
					nOvs++;
					GFFRecord ref = feats.next();
					StrandedFeature.Strand refStrand = ref.getStrand();
					if (refStrand == StrandedFeature.NEGATIVE)
					{
						overlapsRev = true;
					}
					else if (refStrand == StrandedFeature.POSITIVE)
					{
						overlapsFw = true;
					}
					
					String newName = GFFUtils.getGffRecordName(ref);
					if (newName != null) newNames.add(newName);
				}
				
				StrandedFeature.Strand newStrand = StrandedFeature.UNKNOWN;
				if (overlapsRev && !overlapsFw)
				{
					newStrand = StrandedFeature.NEGATIVE;
				}
				else if (overlapsFw && !overlapsRev)
				{
					newStrand = StrandedFeature.POSITIVE;
				}
				target.setStrand(newStrand);
				
				// Transfer name
				if (newStrand != StrandedFeature.UNKNOWN)
				{
					if (newNames.size()>1) System.err.printf("CGI found %d names: %s\n", 
							newNames.size(),ListUtils.excelLine(newNames.toArray(new String[1])));
				}
				
				if (!dontRename && (newNames.size()==1))
				{
					GFFUtils.setGffRecordName(target, newNames.toArray(new String[1])[0]);
				}
				
				
//				Logger.getAnonymousLogger().severe(String.format("GFF %s\t%d overlapping recs\tsawFw=%s\tsawRev=%s\tnewStrand=%s\n", 
//						GFFUtils.gffBetterString(target),nOvs,""+ overlapsFw,"" + overlapsRev, ""+newStrand));
//				System.err.printf("CGI found %d names: %s\n", 
//						newNames.size(),ListUtils.excelLine(newNames.toArray(new String[1])));

				
				outFeats.add_feature(target);
			}
		}
		
		outFeats.write_records(outFn);
		
//		FeatDbParams params = new FeatDbParams();
//		params.addRangeFilter("chr11", 195000, 210000);
//		params.addRangeFilter("chr11", 1195000, 1230000);
//		params.addFeatFilter("knownGene-exon.hg18.clustered0bp.gtf");
//		params.addFeatFilter("DbRepeatMaskerSINE.hg18.chr11.gtf");
//		FeatIterator feats = new FeatIterator(params);
//		while (feats.hasNext())
//		{
//			GFFRecord rec = feats.next();
//			System.out.print(GFFUtils.gffLine(rec));
//		}
		

//		// Setup writer
//		PrintWriter writer = new PrintWriter(new FileOutputStream(String.format("%s.charts.html", outputPrefix)));

	
//			
	}
	
	
}
