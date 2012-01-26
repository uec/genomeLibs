package org.broadinstitute.sting.gatk.uscec.writer;

import java.util.Comparator;
import java.util.PriorityQueue;
import java.util.Set;
import java.util.TreeSet;
import java.util.concurrent.PriorityBlockingQueue;

import org.broad.tribble.util.variantcontext.VariantContext;
import org.broad.tribble.vcf.VCFHeader;
import org.broad.tribble.vcf.VCFWriter;



public abstract class SortingVCFWriterOwnBase implements VCFWriter {

	// The VCFWriter to which to actually write the sorted VCF records
    private VCFWriter innerWriter = null;

    // the current queue of un-emitted records
    private PriorityBlockingQueue<VCFRecord> queue = null;

    // The locus until which we are permitted to write out (inclusive)
    protected Integer mostUpstreamWritableLoc;
    protected static final int BEFORE_MOST_UPSTREAM_LOC = 0; // No real locus index is <= 0

    // The set of chromosomes already passed over and to which it is forbidden to return
    private Set<String> finishedChromosomes = null;

    // Should we call innerWriter.close() in close()
    private boolean takeOwnershipOfInner;

    /**
     * create a local-sorting VCF writer, given an inner VCF writer to write to
     *
     * @param innerWriter        the VCFWriter to write to
     * @param takeOwnershipOfInner Should this Writer close innerWriter when it's done with it
     */
    public SortingVCFWriterOwnBase(VCFWriter innerWriter, boolean takeOwnershipOfInner) {
        this.innerWriter = innerWriter;
        this.queue = new PriorityBlockingQueue<VCFRecord>(1000, new VariantContextComparator());
        this.mostUpstreamWritableLoc = BEFORE_MOST_UPSTREAM_LOC;
        this.finishedChromosomes = new TreeSet<String>();
        this.takeOwnershipOfInner = takeOwnershipOfInner;
    }

    public SortingVCFWriterOwnBase(VCFWriter innerWriter) {
        this(innerWriter, false); // by default, don't own inner
    }

    public void writeHeader(VCFHeader header) {
        innerWriter.writeHeader(header);
    }

    /**
     * attempt to close the VCF file; we need to flush the queue first
     */
    public void close() {
        stopWaitingToSort();

        if (takeOwnershipOfInner)
            innerWriter.close();
    }

    private void stopWaitingToSort() {
        emitRecords(true);
        mostUpstreamWritableLoc = BEFORE_MOST_UPSTREAM_LOC;
    }

    protected void emitSafeRecords() {
        emitRecords(false);
    }

    protected void noteCurrentRecord(VariantContext vc) {
        // did the user break the contract by giving a record too late?
        if (mostUpstreamWritableLoc != null && vc.getStart() < mostUpstreamWritableLoc) // went too far back, since may have already written anything that is <= mostUpstreamWritableLoc
            throw new IllegalArgumentException("Permitted to write any record upstream of position " + mostUpstreamWritableLoc + ", but a record at " + vc.getChr() + ":" + vc.getStart() + " was just added.");
    }

    /**
     * add a record to the file
     *
     * @param vc      the Variant Context object
     * @param refBase the ref base
     */
    public synchronized void add(VariantContext vc, byte refBase) {
        /* Note that the code below does not prevent the successive add()-ing of: (chr1, 10), (chr20, 200), (chr15, 100)
           since there is no implicit ordering of chromosomes:
         */
        VCFRecord firstRec = queue.peek();
        if (firstRec != null && !vc.getChr().equals(firstRec.vc.getChr())) { // if we hit a new contig, flush the queue
            if (finishedChromosomes.contains(vc.getChr()))
                throw new IllegalArgumentException("Added a record at " + vc.getChr() + ":" + vc.getStart() + ", but already finished with chromosome" + vc.getChr());

            finishedChromosomes.add(firstRec.vc.getChr());
            stopWaitingToSort();
        }

        noteCurrentRecord(vc); // possibly overwritten

        queue.add(new VCFRecord(vc, refBase));
        emitSafeRecords();
    }

    private void emitRecords(boolean emitUnsafe) {
        while (!queue.isEmpty()) {
            VCFRecord firstRec = queue.peek();

            // No need to wait, waiting for nothing, or before what we're waiting for:
            if (emitUnsafe || mostUpstreamWritableLoc == null || firstRec.vc.getStart() <= mostUpstreamWritableLoc) {
                queue.poll();
                innerWriter.add(firstRec.vc, firstRec.refBase);
            }
            else {
                break;
            }
        }
    }

    /**
     * Gets a string representation of this object.
     * @return
     */
    @Override
    public String toString() {
        return getClass().getName();
    }

    private static class VariantContextComparator implements Comparator<VCFRecord> {
        public int compare(VCFRecord r1, VCFRecord r2) {
            return r1.vc.getStart() - r2.vc.getStart();
        }
    }

    private static class VCFRecord {
        public VariantContext vc;
        public byte refBase;

        public VCFRecord(VariantContext vc, byte refBase) {
            this.vc = vc;
            this.refBase = refBase;
        }
    }

}
