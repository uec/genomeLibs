package edu.usc.epigenome.genomeLibs;

import java.util.*;
import org.biojava.utils.SmallMap;
import org.biojava.bio.program.gff.GFFEntrySet;
import org.biojava.bio.program.gff.GFFRecord;
import org.biojava.bio.program.gff.SimpleGFFRecord;
import org.biojava.bio.seq.StrandedFeature;
import org.biojava.bio.symbol.*;

public class LocUtils {

	
	public static Location midpoint(Location l)
	{ 
		int mid_i = midpoint_int(l);
		return new PointLocation(mid_i);
	}
	
	
	// Uses rec to get the midpoint
	public static Location midpoint_span(GFFRecord rec, int extent)
	{
		Location loc = LocUtils.record_to_loc(rec);
		return LocUtils.midpoint_span(loc,extent);
	}
	
	public static Location midpoint_span(Location l, int extent)
		{
		int mid_i = midpoint_int(l);
		int s = (mid_i - (int)Math.floor(extent/2));
		s = Math.max(0,s); 
		int e = s + extent - 1;
		return new RangeLocation(s,e);
	}
	
	public static Location expand(GFFRecord rec, int flank)
	{
		return new RangeLocation(rec.getStart() - flank, rec.getEnd() + flank);
	}
	
	public static int midpoint_int(Location l)
	{
		int out = (int)Math.ceil(((double)l.getMax() + (double)l.getMin()) / 2.0);
		return out;
	}
	

	public static int dist_from_mids(Location l1, Location l2)
	{
		return Math.abs(offset_from_mids(l1,l2));
	}

	public static int offset_from_mids(GFFRecord r1, GFFRecord r2)
	{
		return offset_from_mids(record_to_loc(r1), record_to_loc(r2));
	}
	
	public static int offset_from_mids(Location l1, Location l2)
	{
		int m1 = midpoint_int(l1);
		int m2 = midpoint_int(l2);
		
		int dist = m2-m1;
		return dist;
	}

	// If r1 is upstream of r2, result is negative.
	// If r1 is downstream of r2, result is positive.
	// If strand_relative is set to true, we use the directionality of r2
	// instead of chromosomal directionality
	public static int distFromClosestPointsDirectional(GFFRecord r1, GFFRecord r2, boolean strand_relative)
	{
		int abs_dist = dist_from_closest_points(r1,r2);
		int dist = abs_dist;
		
		if (dist != 0) // Not overlapping
		{
			boolean r1_fiveprime = r1.getEnd() < r2.getStart();
			
			// Strand reversal?
			if (strand_relative && (r2.getStrand() == StrandedFeature.NEGATIVE))
			{
				r1_fiveprime = !r1_fiveprime;
			}
			
			dist = (r1_fiveprime) ? -abs_dist : abs_dist;
		}
		
		return dist;
	}
	
	
	public static int dist_from_closest_points(GFFRecord r1, GFFRecord r2)
	{
		return dist_from_closest_points(record_to_loc(r1),record_to_loc(r2));
	}
	

	public static int dist_from_closest_points(Location l1, Location l2)
	{
		return dist_from_closest_points(l1,l2,false);
	}
	
	public static int dist_from_closest_points(Location l1, Location l2, boolean use_shadows)
	{
		if (use_shadows)
		{
			l1 = LocationTools.shadow(l1);
			l2 = LocationTools.shadow(l2);
			return dist_from_closest_points_shadow(l1,l2);
		}
		
		// Iterate through each
		int min_dist = Integer.MAX_VALUE;
		Iterator l1i = l1.blockIterator();
		L1: while (l1i.hasNext())
		{
			Location l1_r = (Location)l1i.next();
			Iterator l2i = l2.blockIterator();
			L2: while (l2i.hasNext())
			{
				Location l2_r = (Location)l2i.next();
				int this_dist = dist_from_closest_points_shadow(l1_r,l2_r);
				min_dist = Math.min(min_dist, this_dist);
				if (min_dist == 0) // the minimum
				{
					break L1;
				}
				
			}

		}
		return min_dist;
	}
	
	public static int dist_from_closest_points_shadow(Location l1, Location l2)
	{
		
		int dist;
		
		if (l1.overlaps(l2))
		{
			dist = 0;
		}
		else
		{
			// Since they are not overlapping, one most be completely
			// ahead of the other.  Just figure out which is ahead,
			// then subtract.
			int max1 = l1.getMax();
			int min2 = l2.getMin();
			int d1 = min2 - max1;
			if (d1 > 0)
			{
				dist = d1;
			}
			else
			{
				int min1 = l1.getMin();
				int max2 = l2.getMax();
				int d2 = min1 - max2;
				dist = d2;
			}
		}	
		
		return dist;
	}
	
	public static Location record_to_loc(GFFRecord rec)
	{
		int s = rec.getStart();
		int e = rec.getEnd();
		Location loc;
		if (s==e)
		{
			loc = new PointLocation(s);
		}
		else
		{
			loc = new RangeLocation(s,e);
		}

		return loc;
	}
	
	// Subtracts the location in loc from the record in rec. Returns
	// a list of new GFFRecord objects.  If the result is contiguous,
	// it returns a single record, otherwise it returns multiple 
	// records representing the disjointed parts
	public static List record_subtract(GFFRecord rec, Location loc)
	{
		List outl = new ArrayList();
		
		
		return outl;
	}
	
	public static List entry_set_loc_list(GFFEntrySet es)
	{
		List outl = new ArrayList(es.size());
		Iterator record_it = es.lineIterator();
		while (record_it.hasNext())
		{
			GFFRecord rec = (GFFRecord)record_it.next();
			Location loc = record_to_loc(rec);
			outl.add(loc);
		}
		
		return outl;
	}
	
	
	
	public static Location entry_set_merged_loc(GFFEntrySet es)
	{
		
		List all_locs = entry_set_loc_list(es);
		Location merge = LocationTools.union(all_locs);
		
		return merge;
	}
	
	public static GFFEntrySet loc_to_entry_set(Location loc)
	{
		GFFEntrySet es = new GFFEntrySet();
		Iterator block_it = loc.blockIterator();
		while (block_it.hasNext())
		{
			Location block = (Location)block_it.next();
			es.add(new SimpleGFFRecord("anon","anon","anon",
					block.getMin(), block.getMax(), 0, null, 0, null, new SmallMap()));
		}
		return es;
	}
	
	// frac = bases_overlapping/denominator
	// denom_type determines the denominator:
	//		-1 denominator is size of r1 (the only one that's not symmetric)
	//		0 denominator is the number of bases covered by both r1 and r2
	//		1 denominator is the size of the smaller between r1 and r2
	//		2 denominator is the size of the larger between r1 and r2
	public static double frac_overlap (GFFRecord r1, GFFRecord r2, int denom_type)
	{
		if (!r1.getSeqName().equalsIgnoreCase(r2.getSeqName())) return 0.0;
		
		Location l1 = LocUtils.record_to_loc(r1);
		Location l2 = LocUtils.record_to_loc(r2);
		
		Location l_int = l1.intersection(l2);
		
		int int_coverage = LocationTools.coverage(l_int);
		
		double denom;
		switch (denom_type) 
		{
		case -1:		denom=LocationTools.coverage(l1); break; 
		case 0:		denom=LocationTools.coverage(l1.union(l2)); break; 
		case 1:		denom=Math.min(LocationTools.coverage(l1),LocationTools.coverage(l2)); break; 
		default:		denom=Math.max(LocationTools.coverage(l1),LocationTools.coverage(l2)); break; 
		}
		
		//System.err.println("\t\t" + l1 + " " + l2 + " = " + int_coverage + " / " + denom);
		
		return (double)int_coverage/denom;
	}
	
	
}
