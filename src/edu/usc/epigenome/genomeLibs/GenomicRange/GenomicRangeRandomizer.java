package edu.usc.epigenome.genomeLibs.GenomicRange;

import java.util.Random;

import edu.usc.epigenome.genomeLibs.GoldAssembly;
import edu.usc.epigenome.genomeLibs.MethylDb.MethylDbUtils;

public class GenomicRangeRandomizer {

	protected String genome = null;
	protected static Random randomGenerator = new Random();
	protected long maxLong = 0;


	public GenomicRangeRandomizer(String inGenome, String inMaskInGtf,
			String inMaskOutGtf, String inMatchDistanceGtf) 
	{
		this.genome = inGenome;

		this.init(inMaskInGtf, inMaskOutGtf, inMatchDistanceGtf);
	}


	protected void init(String inMaskInGtf, String inMaskOutGtf,
			String inMatchDistanceGtf) 
	{
		maxLong = GoldAssembly.getGenomeSize(this.genome);
		
		
		if (inMaskInGtf!=null)
		{
			System.err.println("GenomicRangeRandomizer::maskInGtf not yet implemented");
		}

		if (inMaskOutGtf!=null)
		{
			System.err.println("GenomicRangeRandomizer::maskOutGtf not yet implemented");
		}

		if (inMatchDistanceGtf!=null)
		{
			System.err.println("GenomicRangeRandomizer::matchDistanceGtf not yet implemented");
		}
	}
	
	public GenomicRange generateRandomRange(int length)
	throws Exception
	{
		boolean passes = false;
		GenomicRange gr = null;
		while (!passes)
		{
			// First get our random number
			long globalStartCoord = Math.round((double)this.maxLong * randomGenerator.nextDouble());
			//System.err.println("Generated a random long: " + globalStartCoord);

			gr = GoldAssembly.getLocalCoordinate(globalStartCoord, this.genome);
			int chromLen = GoldAssembly.chromLengthStatic(gr.getChrom(), this.genome);
			gr.setEnd(gr.getStart() + length - 1);
			passes = (gr.getEnd()<chromLen);

			if (!passes)
			{
				System.err.printf("Generate range of length %d (past end of chrom =%s): %s\n", length, passes, gr.toString());
			}

			passes &= MethylDbUtils.CHROMS.contains(gr.getChrom());
			
		}
		return gr;
	}

}
