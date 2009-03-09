package edu.usc.epigenome.genomeLibs;

import java.util.*;
import java.io.*;

import org.biojava.utils.SmallMap;
import org.biojava.bio.seq.StrandedFeature;
import org.biojava.bio.symbol.*;
import org.biojava.bio.program.gff.*;

import org.apache.commons.math.stat.StatUtils;


/*
 * LEGACY CODE FROM PREVIOUS PROJECT.  SOME FUNCTIONS MAY NOT BE RELEVANT
 */

public class GFFUtils {

	
	public static SimpleGFFRecord cloneGffRecord(GFFRecord rec)
	{
		return 
		new SimpleGFFRecord(rec.getSeqName(),rec.getSource(),rec.getFeature(),
				rec.getStart(), rec.getEnd(), rec.getScore(), rec.getStrand(), 
				rec.getFrame(), rec.getComment(), rec.getGroupAttributes());

	}
	
	
	public static GFFRecord tss(GFFRecord rec)
	{
		SimpleGFFRecord tss = new SimpleGFFRecord(rec);
		int tss_coord = (rec.getStrand() == StrandedFeature.NEGATIVE) ? rec.getEnd() : rec.getStart();
		tss.setStart(tss_coord);
		tss.setEnd(tss_coord);
		return tss;
	}
	
	
	public static String gffCsvLineHeader()
	{
		return "name,chrom,score,start,end,source,strand";
	}
	
	public static String gffCsvLine(GFFRecord rec)
	{
		String out = "";
		
		out += GFFUtils.getGffRecordName(rec);
		out += ",";
		out += rec.getSeqName();
		out += ",";
		out += rec.getScore();
		out += ",";
		out += rec.getStart();
		out += ",";
		out += rec.getEnd();
		out += ",";
		out += rec.getSource();
		out += ",";
		StrandedFeature.Strand s = rec.getStrand();
		out += ((s == StrandedFeature.POSITIVE) ? 1 : ((s == StrandedFeature.NEGATIVE) ? -1 : 0));
		
		return out;
	}
	
	public static String gffCsvMinimalLineHeader()
	{
		return "name,chrom,score,start,end,strand";
	}

	public static String gffCsvMinimalLine(GFFRecord rec)
	{
		String out = "";

		out += GFFUtils.getGffRecordName(rec);
		out += ",";
		out += rec.getSeqName();
		out += ",";
		out += rec.getScore();
		out += ",";
		out += rec.getStart();
		out += ",";
		out += rec.getEnd();
		out += ",";
		StrandedFeature.Strand s = rec.getStrand();
		out += ((s == StrandedFeature.POSITIVE) ? 1 : ((s == StrandedFeature.NEGATIVE) ? -1 : 0));

		return out;
	}
	
	public static String htmlTabLineHeader(String session_str)
	{
		
		String out = "";

		if (session_str != null)
		{
			String url = "http://genome.ucsc.edu/cgi-bin/hgSession?";
			url += "hgS_doLoadUrl=submit&hgS_loadUrlName=http://teaview2.hsc.usc.edu/~benb/" + session_str;
			out += "<TR><TD COLSPAN=5><A TARGET=\"USC\" HREF=\"" + url + "\">click for new UCSC browser settings</A></TD></TR>\n";
			out += "<TR><TD COLSPAN=5><BR></TD></TR>\n";
		}
		
		
		out += "<TR>\n";
		
		out += "<TH>name</TH>\n";
		out += "<TH>chrom</TH>\n";
		out += "<TH>start</TH>\n";
		out += "<TH>end</TH>\n";
		out += "<TH>link</TH>\n";
		
		out += "</TR>\n";
		
		
		
		return out;
	}

	public static String htmlTabLine(GFFRecord rec, String tracks_str)
	{
		String out = "";
		
		out += "<TR>\n";
		
		out += "<TD>" + GFFUtils.getGffRecordName(rec) + "</TD>\n";  
		out += "<TD>" + rec.getSeqName() + "</TD>\n";  
		out += "<TD>" + rec.getStart() + "</TD>\n";  
		out += "<TD>" + rec.getEnd() + "</TD>\n";  
		
		out += "<TD>" + recHyperlink(rec, tracks_str) + "</TD>\n";
		
		out += "</TR>\n";
		
		return out;
	}
	
	public static String recHyperlink(GFFRecord rec, String tracks_str)
	{
		String out = "";
		
		int flank = 500;
		
		String url = "http://genome.ucsc.edu/cgi-bin/hgTracks?org=Human&db=hg18&position=";
		url += rec.getSeqName() + ":" + (rec.getStart() - flank) + "-" + (rec.getEnd() + flank);
		if (tracks_str != null)
		{
			url += "&hgt.customText=http://teaview2.hsc.usc.edu/~benb/" + tracks_str;
		}
		
		out = "<A  TARGET=\"USC\" HREF=\"" + url + "\">" + "ucsc browser" + "</A>";

		return out;
	}

	public static String lengthLine(GFFRecord rec)
	{
		String out = "";
		
		out += recLength(rec);
		out += "," + GFFUtils.gffCsvLine(rec);
		
		return out;
	}
	
	public static String scoreLine(GFFRecord rec)
	{
		String out = "";
		
		out += GFFUtils.getGffRecordName(rec);
		out += ",";
		out += rec.getSeqName();
		out += ",";
		out += rec.getStart();
		out += ",";
		out += rec.getEnd();
		out += ",";
		out += rec.getScore();
		
		return out;
	}
	
//	public static String cseLine(GFFRecord rec)
//	{
//		String out = "";
//		
//		out += GFFUtils.gffChromNum(rec);
//		out += ",";
//		out += rec.getStart();
//		out += ",";
//		out += rec.getEnd();
//		
//		return out;
//	}
	
	
	public static int recLength(GFFRecord rec)
	{
		int len = rec.getEnd() - rec.getStart() + 1;
		return len;
	}
	
	public static String gffLine(GFFRecord rec)
	{
		StringWriter strw = new StringWriter(500);
		GFFWriter gffw = new GFFWriter(new PrintWriter(strw));
		gffw.recordLine(rec);
		return strw.toString();
	}
	
	public static String gffBetterString(GFFRecord rec)
	{
		String out = "";
		
		out += GFFUtils.getGffRecordName(rec);
		out += "(";
		out += rec.getSeqName();
		out += ":";
		out += rec.getStart();
		out += "-";
		out += rec.getEnd();
		out += ")";

		return out;
	}
	
	
//	public static int gffChromNum(GFFRecord rec)
//	{
//		int chr = 
//			(new ChromFeatures()).chrom_from_public_str(rec.getSeqName());
//		
//		return chr;
//	}
	
	// -----------------
	// Merge an entry set into a single record
	public static GFFRecord entry_set_merged (GFFEntrySet es, String id_fld, String id)
	{
		return entry_set_merged(es, id_fld, id, new MergeGffParams());
	}
	
	public static GFFRecord entry_set_merged (GFFEntrySet es, String id_fld, String id, MergeGffParams params)
	{
		List<String> l = new ArrayList<String>();
		l.add(id);
		
		Map<String,List<String>> m = (Map<String,List<String>>)new SmallMap();
		m.put(id_fld, l);
		params.f_combine_atts = false;
		return entry_set_merged(es, m, params);
	}
	
	// We concatenate sources if they are different
	public static GFFRecord entry_set_merged (GFFEntrySet es, MergeGffParams params)
	{
		return entry_set_merged(es, null, params);
	}
	
	public static GFFRecord entry_set_merged (GFFEntrySet es, Map<String,List<String>> map, MergeGffParams params)
	{
		Location loc = null;
		try
		{
			loc = LocUtils.entry_set_merged_loc(es);
		}
		catch (Exception e)
		{
			System.err.println("Problem calling entry_set_merged_loc with " + es + "\n" + e.toString());
			return null;
		}
		
		
		
		int num_feats = es.size();
		
		// Get info from the first record
		Iterator<GFFRecord> it = (Iterator<GFFRecord>)es.lineIterator();
		if (!it.hasNext()) return null;
		GFFRecord first_rec = (GFFRecord)it.next();
		
		//System.out.println("Merging" + es.size() + "Records with map: " + map.toString());
		
		// Create new one
		SimpleGFFRecord new_rec =  new SimpleGFFRecord(first_rec);
		
		if (num_feats > 1)
		{
			double[] scores = new double[num_feats];
			scores[0] = first_rec.getScore();
			
			SimpleGFFRecord best_rec = new_rec;
			Set<String> sources = new HashSet<String>();
			sources.add(best_rec.getSource());
			SimpleGFFRecord all_atts_rec = new SimpleGFFRecord(best_rec);
			int on_rec = 1;
			while (it.hasNext())
			{
				SimpleGFFRecord r = (SimpleGFFRecord)it.next();

					if (r.getScore() > best_rec.getScore()) best_rec = r;
					sources.add(r.getSource());

					// Attributes
//					add_all_gffrecord_map_entries(all_atts_rec, r);

					scores[on_rec] = r.getScore();
					
					System.err.println("Including rec with score: " + r.getScore());

					on_rec++;
			}
			new_rec = new SimpleGFFRecord(best_rec);
			
			if (params.f_combine_atts)
			{
				System.err.println("COMBINE_ATTS=" + params.f_combine_atts);
				new_rec.setGroupAttributes(all_atts_rec.getGroupAttributes());
			}
			
			// Put the concatenated sources
			if (params.f_combine_sources)
			{
				new_rec.setSource(ListUtils.excelLine(new ArrayList(sources)));
			}
			
			// The important thing, the start and end
			new_rec.setStart(loc.getMin());
			new_rec.setEnd(loc.getMax());
			
			double score = (params.f_combine_scores) ? StatUtils.mean(scores, 0, on_rec) : on_rec;
			new_rec.setScore(score);
		}

		return new_rec;
	}
	
	public static void addEntries(GFFEntrySet to, GFFEntrySet from)
	{
		Iterator it = from.lineIterator();
		while (it.hasNext())
		{
			to.add((GFFRecord)it.next());
		}
	}
	
	public static double maxScore(GFFEntrySet es)
	{
		Iterator it = es.lineIterator();
		double best = Double.NEGATIVE_INFINITY;
		while (it.hasNext())
		{
			double s = ((GFFRecord)it.next()).getScore();
			if (s>best) best = s;
		}
		return best;
	}
	
//	public static GFFEntrySet combineEntrySets(GFFEntrySet[] ess)
//	{
//		for (int i = 0; i < ess.length; i++)
//		{
//			GFFEntrySet es = 
//		}
//	}
//	
	// ----------
	// Attribute maps

	public static double get_gffrecord_map_entry_double(GFFRecord site, String field)
	{
		String str = (String)get_gffrecord_map_entry_obj(site,field);
		return Double.parseDouble(str);
	}

	public static String get_gffrecord_map_entry(GFFRecord site, String field)
	{
		return (String)get_gffrecord_map_entry_obj(site,field);
	}
	
	public static Object get_gffrecord_map_entry_obj(GFFRecord site, String field)
	{
		Object o = null;
		if (site.getGroupAttributes().containsKey(field))	
		{
			ArrayList al = (ArrayList)site.getGroupAttributes().get(field);
			o = al.get(0);
		}
		return o;
	}
	
	public static List get_gffrecord_map_entries(GFFRecord site, String field)
	{
		List o = null;
		if (site.getGroupAttributes().containsKey(field))	
		{
			ArrayList al = (ArrayList)site.getGroupAttributes().get(field);
			o = al;
		}
		return o;
	}
	
	// Sets if not already set
	public static void checkGffRecordName(GFFRecord site)
	{
		String name = getGffRecordName(site);
		
		if (name == null)
		{
			name = site.getSource() + "-" + site.getSeqName() + "-" + site.getStart();
		}
		
		setGffRecordName(site, name);
	}
	
	
	public static void setGffColor(GFFRecord site, String color)
	{
		put_gffrecord_map_entry(site, "itemRgb", color);
	}

	public static void setGffRecordName(GFFRecord site, String name)
	{
		put_gffrecord_map_entry(site, "gene_id", name);
		put_gffrecord_map_entry(site, "transcript_id", name);
	}
	
	public static String getGffRecordName(GFFRecord site)
	{
		String out = get_gffrecord_map_entry(site, "gene_id");
		return out;
	}

	
	public static void setGffRecordLineNumber(GFFRecord site, String fn, int line_number)
	{
		add_gffrecord_map_entry(site, "ln" + fn, Integer.toString(line_number));
	}
	
	public static int getGffRecordLineNumber(GFFRecord site, String fn)
	{
		String line_str = get_gffrecord_map_entry(site, "ln" + fn);
		int out = Integer.parseInt(line_str);
		
		return out;
	}

	public static List getGffRecordLineNumbers(GFFRecord site, String fn)
	{
		List lines = get_gffrecord_map_entries(site, "ln" + fn);
		return lines;
	}

	public static Set getGffRecordLineNumberIds(GFFRecord rec)
	{
		HashSet s = new HashSet();
		
		Iterator key_it = rec.getGroupAttributes().keySet().iterator();
		while (key_it.hasNext())
		{
			String key = (String)key_it.next();
			if (key.startsWith("ln"))
			{
				String id = key.substring(2,key.length());
				s.add(id);
				// System.err.println("Found key " + key + " - " + id);
			}
		}
		return s;
	}
	
	public static void add_gffrecord_map_entry(GFFRecord site, String field, Object val)
	{
		
		ArrayList al = (ArrayList)site.getGroupAttributes().get(field);
		if (al == null)
		{
			al = new ArrayList();
		}
		al.add(val);
		site.getGroupAttributes().put(field, al);
	}
	
	public static void add_gffrecord_map_entries(GFFRecord site, String field, List vals)
	{
		
		ArrayList al = (ArrayList)site.getGroupAttributes().get(field);
		if (al == null)
		{
			al = new ArrayList();
		}
		al.addAll(vals);
		site.getGroupAttributes().put(field, al);
	}
	
	public static void add_all_gffrecord_map_entries(GFFRecord to, GFFRecord from)
	{
		Map from_map = from.getGroupAttributes();
		Iterator key_it = from_map.keySet().iterator();
		
		while (key_it.hasNext())
		{
			String key = (String)key_it.next();
			ArrayList from_list = (ArrayList)from_map.get(key);
			add_gffrecord_map_entries(to, key, from_list);
		}
	}

	public static void put_gffrecord_map_entry(GFFRecord site, String field, Object val)
	{
		ArrayList al = new ArrayList();
		al.add(val);
		site.getGroupAttributes().put(field, al);
	}
	
	public static boolean conservedAttribute(GFFRecord rec)
	{
		String cons_str = GFFUtils.get_gffrecord_map_entry(rec, "conserved");
		
		boolean conserved = false;
		if (cons_str != null)
		{
			conserved = cons_str.equals("+");
		}

		return conserved;
	}
	

	// -----------
	// Location processing
	

	
	// -----------
	// Sorting, diffing
	
	public static GFFEntrySet sortEntrySet(GFFEntrySet in)
	{
		return sortEntrySet(in, GFFRecord.NATURAL_ORDER);
	}
	
	public static GFFEntrySet sortEntrySet(GFFEntrySet in, Comparator comp)
	{

		// Load ordered list
		List list = new ArrayList(in.size()); 
		Iterator in_it = in.lineIterator();
		while (in_it.hasNext())
		{
			list.add(in_it.next());
		}

		// Sort
		Collections.sort(list,comp);
		
		// Now load the new entry set
		GFFEntrySet out = new GFFEntrySet();
		Iterator list_it = list.iterator();
		while (list_it.hasNext())
		{
			out.add((GFFRecord)list_it.next());
		}
		
		return out;
	}
	
	// Removes all entries in B from A
	//
	// If A=1,2,3,4
	// and B=2,4
	//
	// returns C=1,3
	public static GFFEntrySet subtractSet(GFFEntrySet a, GFFEntrySet b)
	{
		return a.filter( new InEntrySetFilter(b, true) );
	}
	
	/***** 
	 * Cluster finding
	 * 
	 */
	public static GFFEntrySet simpleClusters(GFFEntrySet in_es, int max_flank, 
			MergeGffParams params)
	{
		return simpleClusters(in_es, max_flank, 0, params);
	}
	
	public static GFFEntrySet simpleClusters(GFFEntrySet in_es, int max_flank, 
			int max_clust_size, MergeGffParams params)
	{
		int num_in = in_es.size();
		if (num_in==0) return new GFFEntrySet();
		if (max_clust_size <= 0) max_clust_size = Integer.MAX_VALUE;
		
		//System.err.println("simpleClusters("+num_in+","+max_flank+","+min_features+")");
		
		// Make sure records are sorted
		in_es = sortEntrySet(in_es);
		
		// Output list
		GFFEntrySet clusts = new GFFEntrySet();
		
	
		// Loop through records one by one
		GFFEntrySet clust = new GFFEntrySet();
		Iterator in_it = in_es.lineIterator();
		int clust_end = Integer.MIN_VALUE;
		int clust_start = Integer.MAX_VALUE;
		int rec_num = 0;
		while (in_it.hasNext())
		{
			rec_num++;
			GFFRecord rec = (GFFRecord)in_it.next();
			
			// Decide if this record should belong to the current cluster
			if (includeInSimpleCluster(clust, rec, in_es, max_flank, max_clust_size, clust_start, clust_end, params))
			{
				clust.add(rec);
				clust_end = Math.max(clust_end, rec.getEnd());
				clust_start = Math.min(clust_start, rec.getStart());
			}
			else
			{
				//System.err.println("\tStarting a new cluster");
				if (clust.size()>=params.f_min_features)
				{
					clusts.add(entry_set_merged(clust, params));
				}
				
				// Start a new cluster
				clust = new GFFEntrySet();
				clust_end = Integer.MIN_VALUE;
				clust_start = Integer.MAX_VALUE;
			
				
				if (rec.getScore() >= params.f_min_score)
				{
					clust.add(rec);
					clust_start = rec.getStart();
					clust_end = rec.getEnd();
				}
			}

			if ((rec_num % 5000) == 0)
			{
				System.err.println("Working on rec " + rec_num + "/" + num_in + "\t" + LocUtils.record_to_loc(rec));
			}
		}
		
		// Finish off the last cluster
		if (clust.size()>=params.f_min_features)
		{
			clusts.add(entry_set_merged(clust, params));
		}

		return clusts;
	}

	protected static boolean includeInSimpleCluster(GFFEntrySet clust, GFFRecord rec,
			GFFEntrySet all_records, int max_flank, int max_clust_size, int clust_start, int clust_end,
			MergeGffParams params)
	{
		if (rec.getScore() < params.f_min_score) return false;
		
		
		// If the cluster hasn't started yet, add it.
		if (clust.size() == 0) return true;

		int rec_start = rec.getStart();
		int rec_end = rec.getEnd();
		
		// If the current cluster is too big, end it.
		if ( (rec_end - clust_start) > max_clust_size ) return false;
		
		// Location es_loc = LocUtils.entry_set_merged_loc(clust); // Slow
		// int dist = LocUtils.dist_from_closest_points(es_loc, rec_loc, true);

		int dist = rec_start - clust_end;
		//Location rec_loc = LocUtils.record_to_loc(rec);
		
		// System.err.println("\tClust size=" + clust.size() + "\tdist=" + dist + "\tmax_flank=" + max_flank);
		return (dist <= max_flank);
	}
		
}
