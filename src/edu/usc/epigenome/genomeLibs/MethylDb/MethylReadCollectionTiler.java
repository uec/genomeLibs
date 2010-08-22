package edu.usc.epigenome.genomeLibs.MethylDb;

import java.io.PrintWriter;
import java.util.Iterator;
import java.util.Set;
import java.util.TreeSet;

public class MethylReadCollectionTiler extends MethylReadCollection {

	protected Set<Integer> cpgPositions = null;
	
	
	public MethylReadCollectionTiler(MethylDbQuerier inParams) throws Exception {
		super(inParams);
		init();
	}

	/* (non-Javadoc)
	 * @see edu.usc.epigenome.genomeLibs.MethylDb.MethylReadCollection#init()
	 */
	@Override
	protected void init() {
		super.init();
		cpgPositions = new TreeSet<Integer>();
	}
	
	
	/******** INFO ************/



	/******** POPULATING ************/

	/* (non-Javadoc)
	 * @see edu.usc.epigenome.genomeLibs.MethylDb.MethylReadCollection#addCpg(edu.usc.epigenome.genomeLibs.MethylDb.Cpg)
	 */
	@Override
	public void addCpg(Cpg cpg) {
		super.addCpg(cpg);
		cpgPositions.add(cpg.chromPos);
	}

	
	/******** OUTPUT ************/
	
	public void writeTiling(PrintWriter pw)
	{
		for (Integer i : cpgPositions)
		{
			pw.printf("Pos %d\n", i);
		}
		
		// Make a sorted list of reads
		Set<MethylRead> reads = this.getSortedReads();
		for (MethylRead read : reads)
		{
			pw.printf("Read: %s\n",read.toString());
		}
		
	}

	
	
	
}
