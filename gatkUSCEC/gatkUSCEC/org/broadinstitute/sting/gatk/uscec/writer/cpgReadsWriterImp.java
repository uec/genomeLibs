package org.broadinstitute.sting.gatk.uscec.writer;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;

public class cpgReadsWriterImp extends FormatWriterBase {

	public cpgReadsWriterImp(File location) {
		super(location);
		// TODO Auto-generated constructor stub
	}

	public cpgReadsWriterImp(OutputStream output) {
		super(output);
		// TODO Auto-generated constructor stub
	}

	public cpgReadsWriterImp(File location, OutputStream output) {
		super(location, output);
		// TODO Auto-generated constructor stub
	}

	@Override
	public void add(genomeObject obj) {
		// TODO Auto-generated method stub
		String readsLine = String.format("%s\t%d\t%c\t%d\t%c\t%s\t%d\n",obj.getChr(), obj.getStart(), ((cpgReads) obj).getMethyStatus(), ((cpgReads) obj).getbaseQ(), ((cpgReads) obj).getstrand(), ((cpgReads) obj).getReadID(), ((cpgReads) obj).getEncryptID());
    	try {
			mWriter.write(readsLine);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	@Override
	public void addHeader(Object o) {
		// TODO Auto-generated method stub

	}

}
