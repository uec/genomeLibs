package org.broadinstitute.sting.gatk.uscec.writer;

public class cpgReads implements genomeObject {

	private String chr;
	private int genomeLoc;
	private char methyStatus;
	private byte baseQ;
	private char strand;
	private String readID;
	
	
	public cpgReads(String chr, int genomeLoc, char methyStatus, byte baseQ, char strand, String readID){
		this.chr = chr;
		this.genomeLoc = genomeLoc;
		this.methyStatus = methyStatus;
		this.baseQ = baseQ;
		this.strand = strand;
		this.readID = readID;
	}
	
	
	public char getMethyStatus(){
		return this.methyStatus;
	}
	
	public byte getbaseQ(){
		return this.baseQ;
	}
	
	public char getstrand(){
		return this.strand;
	}
	
	public String getReadID(){
		return this.readID;
	}
	
	@Override
	public int getStart() {
		// TODO Auto-generated method stub
		return this.genomeLoc;
	}

	@Override
	public String getChr() {
		// TODO Auto-generated method stub
		return this.chr;
	}

}
