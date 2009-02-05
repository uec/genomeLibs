package edu.usc.epigenome.genomeLibs.AlignmentPos;

import edu.usc.epigenome.genomeLibs.ReadPos.ReadPosOptions;

public class AlignmentPosOptions extends ReadPosOptions {

	public String f_genome = "hg18";
	public int maxIdentical = 0;
	public boolean onlyFirstCycle = false; 
	public boolean trackSnps = true; // Will force iterators to use AlignmentPosSnps objects
	public boolean trackBisulfiteConversion = false; // Will force iterators to use AlignmentPosBisulfiteConvSnps objects

}
