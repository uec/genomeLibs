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
				String name = String.format("%s-%d-%d-cpg%d",
						walkers[0].getCurChr(), CpgWalker.windStart(wind), CpgWalker.windEnd(wind), pos);
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

		List<String> fields = new ArrayList<String>(20);
		
		int pos = cpgs[0].chromPos;
		
		fields.add(this.getCurChr());
		fields.add(Integer.toString(pos));
		fields.add(Integer.toString(pos));
		fields.add(name);
		fields.add(Character.toString(cpgs[0].getStrand().getToken()));
		
		for (double meth : meths)
		{
			fields.add(Double.toString(meth));
		}
		
		ListUtils.setDelim("\t");
		pw.println( ListUtils.excelLine(fields) );
		
	}

}
