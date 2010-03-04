package edu.usc.epigenome.genomeLibs;

public interface TabularOutput {

	public String toCsvStr() throws Exception;
	public String headerStr() throws Exception;
	
}
