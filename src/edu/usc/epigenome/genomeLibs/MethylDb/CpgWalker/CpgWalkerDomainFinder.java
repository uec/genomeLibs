package edu.usc.epigenome.genomeLibs.MethylDb.CpgWalker;

import java.util.ArrayList;
import java.util.List;

import edu.usc.epigenome.genomeLibs.GenomicRange.*;

public class CpgWalkerDomainFinder extends CpgWalker {

	List<GenomicRange> domains = new ArrayList<GenomicRange>();
	
	public CpgWalkerDomainFinder(CpgWalkerParams inWalkParams) {
		super(inWalkParams);
		// TODO Auto-generated constructor stub
	}

	@Override
	protected void processWindow() {
		// TODO Auto-generated method stub

	}

}
