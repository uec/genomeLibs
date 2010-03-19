package edu.usc.epigenome.genomeLibs.MethylDb.CpgWalker;

import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;

import edu.usc.epigenome.genomeLibs.ListUtils;
import edu.usc.epigenome.genomeLibs.MethylDb.Cpg;

public class CpgWalkerMultisampleMethRange extends CpgWalkerMultisample {

	PrintWriter pw = null;
	double minMeth = Double.NaN;
	double maxMeth = Double.NaN;
	
	int lastPosOutput = -1;
	int curStartPos = -1;
	int nextWindId = 1;
	
	public CpgWalkerMultisampleMethRange(CpgWalkerParams inWalkParams,
			int numSamples, PrintWriter inPw, double inMinMeth, double inMaxMeth) {
		super(inWalkParams, numSamples);

		pw = inPw;
		minMeth = inMinMeth;
		maxMeth = inMaxMeth;
		
	}

	
	
	@Override
	public void reset() {
		super.reset();
		lastPosOutput = -1;
		curStartPos = -1;
		System.err.println("CpgWalker multisample reset()");
	}



	/* (non-Javadoc)
	 * @see edu.usc.epigenome.genomeLibs.MethylDb.CpgWalker.CpgWalkerMultisample#processWindow()
	 * 
	 * This has a specialized purpose for now.  The point is to identify all CpGs lying in windows
	 * where one (or more) of the samples fall within the meth range, and one (or more) of the
	 * samples fall outside of it.  Then we output all CpGs in these windows. We also output
	 * a flag if the infividual CpG itself follows the trend.
	 */
	@Override
	public void processWindow() {
		
		boolean[] inside = new boolean[this.numSamples()];
		boolean[] outside = new boolean[this.numSamples()];
		int totalInside = 0;
		int totalOutside = 0;
		for (int i = 0; i < this.numSamples(); i++)
		{
			CpgWalker walker = walkers[i];
			double methMin = walker.methSummarizer.getValsMin();
			double methMax = walker.methSummarizer.getValsMax();
			inside[i] = (methMin>=this.minMeth) && (methMax<=this.maxMeth);
			outside[i] = (methMax<this.minMeth) || (methMin>this.maxMeth);
			if (inside[i]) totalInside++;
			if (outside[i]) totalOutside++;
		}
		
//		boolean varies = (totalInside>0) && (totalInside<numSamples());
		
		// This is the very strict definition
		boolean varies = (totalInside>0) && (totalOutside>0);
		
//		if (totalInside>0)
//		{
//			System.err.printf("Wind: %s\tnumInside=%d\tnumOutside=%d\tvaries=%s\n",
//					this.lastProcessedWindStr(false),totalInside,totalOutside,varies);
//		}
		
		if (varies) outputWindow();
	}



	protected void outputWindow() {

		double methDiff = walkers[0].methSummarizer.getValMean(false) -  walkers[1].methSummarizer.getValMean(false);// This is not a good measure for more than 2 samples
		
		
		// All of the processed windows should have identical CpGs
		Object[] winds = new Object[numSamples()];
		for (int i = 0; i < numSamples(); i++)
		{
			winds[i] = walkers[i].getLastProcessedWindow();
		}
		
		List<Cpg> wind1 = (List<Cpg>)winds[0];
		int nCpgs = wind1.size();
		
		for (int i = 0; i < nCpgs; i++)
		{
			int pos = wind1.get(i).chromPos;
			
			if (pos>lastPosOutput)
			{
				int totalInside = 0;
				Cpg[] cpgs = new Cpg[numSamples()];
				double[] meths = new double[numSamples()]; 
				for (int s = 0; s < numSamples(); s++)
				{
					Cpg cpg = ((List<Cpg>)winds[s]).get(i);
					int posi = cpg.chromPos;
					double methi = cpg.fracMeth(true);
					boolean inside = (methi>=this.minMeth) && (methi<=this.maxMeth);
					if (inside) totalInside++;

					if (posi != pos)
					{
						System.err.printf("\tpositions do not match %d != %d\n",posi,pos);
						System.exit(1);
					}
					//				else
					//				{
					//					System.err.printf("\tpositions DO match %d == %d\n",posi,pos);
					//				}
					
					cpgs[s] = cpg;
					meths[s] = methi;
				}
				
				
				List<Cpg> wind = walkers[0].getLastProcessedWindow();
				int ws = CpgWalker.windStart(wind);
				int we = CpgWalker.windEnd(wind);
				
				// If the start pos of this window is before the last pos output,
				// we merge this window with the prior one.
				if (ws>lastPosOutput) curStartPos = ws;
				
				String name = String.format("%s-s%d",
						walkers[0].getCurChr(), curStartPos);
//						((methDiff>=0.0)?"+":""),methDiff
				boolean varies = (totalInside>0) && (totalInside<numSamples());
				varies = true;
				if (varies)
				{
					outputCpg(cpgs, meths,name);
				}

				// And update last pos seen
				lastPosOutput = pos;

			}

		
		}
				
		
	}



	private void outputCpg(Cpg[] cpgs, double[] meths, String name) {

//		chr11   DEFAULT exon    259328  259328  14      +       .       gene_id "chr11-s259328-d+0.61"; transcript_id "chr11-s259328-d+0.61";		
		
		//!!!!! WE HAVE TO MAKE A SPECIAL CASE FOR LISTER, FOR WHICH SOME
		// HAVE A FAKE 1 READ.
		boolean suffReads = true;
		int minCtReads = Integer.MAX_VALUE; 
		StringBuffer methSec = new StringBuffer(500);
		for (int i = 0; i < cpgs.length; i++)
		{
			Cpg cpg = cpgs[i];
			if (cpg.totalReads < this.walkParams.minReadsForOutput) suffReads = false;
			minCtReads = Math.min(minCtReads, cpg.cReads+cpg.tReads);
			methSec.append(String.format("meth%d \"%.2f\"; ", i, meths[i]));
		}
		if (!suffReads) return;
		
		List<String> fields = new ArrayList<String>(20);
		
		int pos = cpgs[0].chromPos;
		
		fields.add(this.getCurChr());
		fields.add("USC");
		fields.add("exon");
		fields.add(Integer.toString(pos));
		fields.add(Integer.toString(pos));
		fields.add(Integer.toString(minCtReads)); // score
		fields.add(Character.toString(cpgs[0].getStrand().getToken()));
		fields.add(".");
		fields.add(String.format("gene_id \"%s\"; transcript_id \"%s\"; %s",name,name,methSec.toString()));
		
		ListUtils.setDelim("\t");
		pw.println( ListUtils.excelLine(fields) );
		
	}

}
