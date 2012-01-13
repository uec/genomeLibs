package org.broadinstitute.sting.gatk.uscec.bisulfitesnpmodel;


import org.broadinstitute.sting.gatk.uscec.writer.SortingVCFWriterOwn;

public class SortingTcgaVCFWriter extends SortingVCFWriterOwn {

	protected TcgaVCFWriter tcgaInnerWriter = null;
	public SortingTcgaVCFWriter(TcgaVCFWriter innerWriter,
			int maxCachingStartDistance, boolean takeOwnershipOfInner) {
		super(innerWriter, maxCachingStartDistance, takeOwnershipOfInner);
		tcgaInnerWriter = innerWriter;
		// TODO Auto-generated constructor stub
	}

	public SortingTcgaVCFWriter(TcgaVCFWriter innerWriter,
			int maxCachingStartDistance) {
		super(innerWriter, maxCachingStartDistance);
		tcgaInnerWriter = innerWriter;
		// TODO Auto-generated constructor stub
	}

	public TcgaVCFWriter getInnerWriter(){
		return this.tcgaInnerWriter;
	}
	

	
}
