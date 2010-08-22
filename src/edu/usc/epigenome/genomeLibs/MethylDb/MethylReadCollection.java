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
		//System.err.printf("Adding Cpg pos=%d, numReads=%d, meth=%.2f\n%s\n", cpg.chromPos, cpg.totalReads, cpg.fracMeth(true),cpg.toReadFormatString());
		
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
					//System.err.printf("Pos %d, Read exists\n",cpg.chromPos);
				}

				boolean meth = cpgRead.meth(this.params.useNonconversionFilter);
				read.addPosition(cpg.chromPos, meth);
//				System.err.printf("Read %d has %d CpGs (meth=%.2f)\n", cpgRead.readId, read.numTotal(), read.fracMeth());
			}
		}
	}
	
	public Set<MethylRead> getSortedReads()
	{
		TreeSet<MethylRead> out = new TreeSet<MethylRead>(readMapByReadId.values());
		return out;
	}

	public int numReads() {
		return readMapByReadId.size();
	}
}
