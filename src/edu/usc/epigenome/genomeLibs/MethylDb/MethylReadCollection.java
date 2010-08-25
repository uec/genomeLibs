package edu.usc.epigenome.genomeLibs.MethylDb;

import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;

public class MethylReadCollection {
	
	protected Map<Integer,MethylRead> readMapByReadId = null;
	protected MethylDbQuerier params = null;


	public MethylReadCollection(MethylDbQuerier inParams) throws Exception 
	{
		params = inParams;
		this.init();
	}
	
	protected void init()
	{
		readMapByReadId = new HashMap<Integer,MethylRead>();
	}

	
	public void populate(Iterator<Cpg> cpgIt)
	{
		Cpg cpg = null;
		while (cpgIt.hasNext())
		{
			cpg = cpgIt.next();
			this.addCpg(cpg);
		}
	}
	
	public void addCpg(Cpg cpg)
	{
		addCpg(cpg,false);
	}
	
	/**
	 * @param cpg
	 * @param collapseRevStrandPositions If true, we subtract one from the position of all reverse strand Cpgs, making them the same as their FW strand counterparts 
	 * 
	 */
	public void addCpg(Cpg cpg, boolean collapseRevStrandPositions)
	{
		//System.err.printf("Adding Cpg pos=%d, numReads=%d, meth=%.2f\n%s\n", cpg.chromPos, cpg.totalReads, cpg.fracMeth(true),cpg.toReadFormatString());
		
		int newChromPos = (collapseRevStrandPositions && cpg.negStrand) ? cpg.chromPos-1 : cpg.chromPos;
		cpg.chromPos = newChromPos;
		
		Iterator<CpgRead> readIt = cpg.getReads().values().iterator(); 
		while (readIt.hasNext())
		{
			CpgRead cpgRead = readIt.next();
			if (cpgRead.validCg(this.params.useNonconversionFilter))
			{
				Integer key = new Integer(cpgRead.readId);
				MethylRead read = readMapByReadId.get(key);
				if (read==null)
				{
					//System.err.printf("Creating read for read %d\n",cpgRead.readId);
					read = new MethylRead(cpg.getStrand(), cpgRead.readId);
					readMapByReadId.put(key, read);
				}
				else
				{
					//System.err.printf("Pos %d, Read exists\n",newChromPos);
				}

				
				boolean meth = cpgRead.meth(this.params.useNonconversionFilter);
				read.addPosition(newChromPos, meth);
//				System.err.printf("Read %d has %d CpGs (meth=%.2f)\n", cpgRead.readId, read.numTotal(), read.fracMeth());
			}
		}
	}
	
	public TreeSet<MethylRead> getSortedReads()
	{
		TreeSet<MethylRead> out = new TreeSet<MethylRead>(readMapByReadId.values());
		return out;
	}

	public int numReads() {
		return readMapByReadId.size();
	}
}
