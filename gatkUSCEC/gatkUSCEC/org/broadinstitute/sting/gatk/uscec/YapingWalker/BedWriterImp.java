package org.broadinstitute.sting.gatk.uscec.YapingWalker;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;

import org.broad.tribble.util.variantcontext.VariantContext;
import org.broadinstitute.sting.utils.GenomeLoc;

public class BedWriterImp extends FormatWriter {

	public BedWriterImp(File location) {
		super(location);
		// TODO Auto-generated constructor stub
	}

	public BedWriterImp(OutputStream output) {
		super(output);
		// TODO Auto-generated constructor stub
	}

	public BedWriterImp(File location, OutputStream output) {
		super(location, output);
		// TODO Auto-generated constructor stub
	}

	@Override
	public void add(GenomeLoc loc, double value) {
    	String bedLine = String.format("%s\t%d\t%d\t%.2f\n",loc.getContig(), loc.getStart(), loc.getStop(), value);
    	try {
			mWriter.write(bedLine);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
    }
    
	@Override
    public void add(String contig, long start, long end, double value) {
    	String bedLine = String.format("%s\t%d\t%d\t%.2f\n",contig, start, end, value);
    	try {
			mWriter.write(bedLine);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
    }

	@Override
	public void addHeader() {
		// TODO Auto-generated method stub

	}
	
	
	private static class bedRecord {
        public VariantContext vc;
        public byte refBase;

        public bedRecord(VariantContext vc, byte refBase) {
            this.vc = vc;
            this.refBase = refBase;
        }
    }

}
