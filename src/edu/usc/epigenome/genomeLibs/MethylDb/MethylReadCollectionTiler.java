package edu.usc.epigenome.genomeLibs.MethylDb;

import java.io.PrintWriter;
import java.util.Arrays;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;

import javax.xml.transform.OutputKeys;
import javax.xml.transform.Transformer;
import javax.xml.transform.TransformerException;
import javax.xml.transform.TransformerFactory;
import javax.xml.transform.dom.DOMSource;
import javax.xml.transform.stream.StreamResult;

import org.apache.batik.dom.svg.SVGDOMImplementation;
import org.biojava.bio.seq.StrandedFeature;
import org.w3c.dom.DOMImplementation;
import org.w3c.dom.Document;
import org.w3c.dom.Element;

import com.sun.tools.javac.code.Attribute.Array;

import edu.usc.epigenome.genomeLibs.MiscUtils;

/**
 * @author benb
 *
 */
public class MethylReadCollectionTiler extends MethylReadCollection {

	public final static int MAX_LEVELS = 1000;
	
	protected TreeSet<Integer> cpgPositions = null;
	
	public int imageWidth = 400;
	public int imageHeight = 400;
	
	public int rowHeight = 20;
	public int boxWidth = 3;
	public float boxHeightRelative = 0.8f;
	public float lineHeightRelative = 0.1f;
	public double BOX_BUFFER = 3.0;
	
	public String methColor = "red";
	public String unmethColor = "green";
	
	public String fwStrandLineColor = "black";
	public String revStrandLineColor = "blue";
	
	boolean remapCollisions = true;
	boolean collapseRevStrandPositions = true;
	
	public double[] methTiers = { Double.NEGATIVE_INFINITY, 0.001, 0.999,Double.POSITIVE_INFINITY };
	
	
	TreeSet<MethylRead>[] levels = null;
	
	
	
	public MethylReadCollectionTiler(MethylDbQuerier inParams) throws Exception {
		super(inParams);
		init();
	}

	/* (non-Javadoc)
	 * @see edu.usc.epigenome.genomeLibs.MethylDb.MethylReadCollection#init()
	 */
	@Override
	protected void init() {
		super.init();
		cpgPositions = new TreeSet<Integer>();
		this.clearLevels();
	}

	protected void clearLevels()
	{
		levels = new TreeSet[MAX_LEVELS];
	}
	
	/******** INFO ************/

	public int numLevels()
	{
		int nL = 0;
		for (TreeSet<MethylRead> t : this.levels)
		{
			if (t == null) return nL;
			nL++;
		}
		return nL;
	}
	

	/******** POPULATING ************/

	/* (non-Javadoc)
	 * @see edu.usc.epigenome.genomeLibs.MethylDb.MethylReadCollection#addCpg(edu.usc.epigenome.genomeLibs.MethylDb.Cpg)
	 */
	@Override
	public void addCpg(Cpg cpg) {
		super.addCpg(cpg, this.collapseRevStrandPositions);
		cpgPositions.add(cpg.chromPos);
	}

	
	/* (non-Javadoc)
	 * @see edu.usc.epigenome.genomeLibs.MethylDb.MethylReadCollection#addCpg(edu.usc.epigenome.genomeLibs.MethylDb.Cpg, boolean)
	 */
	@Override
	public void addCpg(Cpg cpg, boolean collapseRevStrandPositions) {
		super.addCpg(cpg, collapseRevStrandPositions);
		cpgPositions.add(cpg.chromPos);
	}

	/******** OUTPUT ************/
	
	public Map<Integer,Double> cpgPositionMapping()
	{
		TreeMap<Integer,Double> map = new TreeMap<Integer,Double>();
		
		// halfWidth puts them right next to each other
		double buffer = BOX_BUFFER + (double)Math.round((float)this.boxWidth/2.0f);
		
		int onCpg = 0;
		for (Integer i : this.cpgPositions)
		{
			
			Double lastPos = (onCpg==0) ? Integer.MIN_VALUE : map.get(map.lastKey());
			double curSpace = (double)i - lastPos;
			if (curSpace<buffer)
			{
				double moveBy = buffer-curSpace;
				System.err.printf("%d less than %.2f from %.2f, moving first %d entries by %.2f (Cpg %d)\n",
						i,buffer,lastPos,map.size(),moveBy,onCpg);
				for (Integer j : map.keySet())
				{
					map.put(j, new Double(map.get(j)-moveBy));
				}
			}
			else
			{
				System.err.printf("%d greater than %.2f from %.2f, not moving\n",
						i,buffer,lastPos);
				
			}
			
			map.put(i, new Double(i));
			onCpg++;
		}
		
		debugTree(map);
		
		
		return map;
	}
	
	public static void debugTree(Map<Integer,Double> map)
	{
		System.err.println("Map contents:");
		for (Integer key : map.keySet())
		{
			Double val = map.get(key);
			System.err.printf("%d (%s) --> (%.2f) %s\n", key, key.getClass(), val, val.getClass());
		}
	}
	
	public void writeTiling(PrintWriter pw)
	{
//		for (Integer i : cpgPositions)
//		{
//			pw.printf("Pos %d\n", i);
//		}
		
		// Make a sorted list of reads
		TreeSet<MethylRead> reads = this.getSortedReads();
		
		// Get the position map
		Map<Integer,Double> positionMap = (this.remapCollisions) ? this.cpgPositionMapping() : null;

		
		// Start the document
		DOMImplementation impl = SVGDOMImplementation.getDOMImplementation();
		String svgNS = SVGDOMImplementation.SVG_NAMESPACE_URI;
		Document doc = impl.createDocument(svgNS, "svg", null);

		// Get the root element (the 'svg' element).
		Element svgRoot = doc.getDocumentElement();

		// Make the level structure
		int startLevel = 0;
		int nLevels = this.methTiers.length - 1;
		if (nLevels <= 1)
		{
			System.err.printf("Forget to add a final element to MethylReadCollectionTiler.methTiers, adding +Infinity\n");
			this.methTiers = new double[2];
			this.methTiers[0] = Double.NEGATIVE_INFINITY;
			this.methTiers[1] = Double.POSITIVE_INFINITY;
		}
		for (int i = 0; i<nLevels; i++)
		{
			double tierLow = methTiers[i];
			double tierHigh = methTiers[i+1];
		
			this.clearLevels();
			for (MethylRead read : reads)
			{
				//			pw.printf("Read: %s\n",read.toString());

				double meth = read.fracMeth();
				boolean add = (meth>tierLow) && (meth<=tierHigh);
				if (add) addRead(read);
			}

			boolean upsideDown = (nLevels>1) && (i==0);
			int newLevels = levelsToDoc(doc, svgRoot, svgNS, startLevel, cpgPositions.first().intValue(), cpgPositions.last().intValue(), upsideDown, positionMap);
			startLevel += newLevels;
		}
		// And write the XML
		writeDoc(doc,pw);
	
		// Clean up
		levels = null;		
	}
		
	
	/**
	 * Expects this.levels to contain reads.
	 * @param doc
	 * @param svgRoot
	 * @param svgNS
	 * @param levelOffset
	 * @param firstCoord
	 * @param lastCoord
	 * @param upsideDown
	 * @param cpgPositionMap If not null, we transform cpg positions based on this map (oldPosition->newPosition)
	 * @return The number of levels used
	 */
	public int levelsToDoc(Document doc, Element svgRoot, String svgNS, int levelOffset, int firstCoord, int lastCoord, boolean upsideDown)
	{
		return levelsToDoc(doc,svgRoot,svgNS,levelOffset,firstCoord,lastCoord,upsideDown,null);	
	}
	
	public int levelsToDoc(Document doc, Element svgRoot, String svgNS, int levelOffset, int firstCoord, int lastCoord, boolean upsideDown, Map<Integer,Double> cpgPositionMap)
	{
		
		TreeSet<MethylRead>[] curLevels = Arrays.copyOf(this.levels,numLevels());
		if (upsideDown) MiscUtils.reverseArray(curLevels);
		

		// Set the width and height attributes on the root 'svg' element.
		svgRoot.setAttributeNS(null, "width", "400");
		svgRoot.setAttributeNS(null, "height", "450");

		int rowHeightHalf = Math.round((float)this.rowHeight / (float)2.0);
		int boxHeight = Math.round(this.boxHeightRelative * (float)this.rowHeight);
		int boxHeightHalf = Math.round(0.5f * this.boxHeightRelative * (float)this.rowHeight);
		int lineHeight = Math.round(this.lineHeightRelative * (float)this.rowHeight);
		int lineHeightHalf = Math.round(0.5f * this.lineHeightRelative * (float)this.rowHeight);
		int boxWidthHalf = (int)Math.floor(0.5f * (float)this.boxWidth);
				
		boolean finished = false;
		int out = curLevels.length;
		for (int l = 0; (l<curLevels.length)&(!finished); l++)
		{
			if (curLevels[l] == null)
			{
				finished = true;
				out = l;
			}
			else
			{
				
				
//				public int rowHeight = 20;
//				public double boxHeightRelativ = 0.8;
//				public double lineHeightRelative = 0.1;
				
				
				// Get the Y information
				int rowCenter = ((l+levelOffset) * this.rowHeight) + rowHeightHalf;
				int lineYStart = rowCenter - lineHeightHalf;
				int lineYEnd = rowCenter + lineHeightHalf;
				
				int boxYStart = rowCenter - boxHeightHalf;
//				int boxYEnd = rowCenter + boxHeightHalf;
				
				
				
				System.err.printf("Level %d: ", l);
				Iterator<MethylRead> rIt = curLevels[l].iterator();
				while (rIt.hasNext())
				{
					MethylRead r = rIt.next();
					double readStartPos = r.startPos();
					double readEndPos = r.endPos();
					if (cpgPositionMap != null)
					{
						Double obj;
						
						obj = cpgPositionMap.get(new Integer((int)readStartPos));
						if (obj == null) System.err.printf("Can't find read start position %.2f in cpgPosition map\n",readStartPos);
						readStartPos = obj.intValue();

						obj = cpgPositionMap.get(new Integer((int)readEndPos));
						if (obj == null) System.err.printf("Can't find read end position %.2f in cpgPosition map\n",readEndPos);
						readEndPos = obj.intValue();
					}
					
					System.err.printf("%.1f-%.1f (%d,%d), ",(readStartPos-firstCoord),(readEndPos-firstCoord),r.numTotal(),(int)(100.0*r.fracMeth()));
					
					// Get the X information
					double lineXStart = readStartPos - (double)firstCoord; 
					double lineXEnd = readEndPos - (double)firstCoord; 
					double lineWidth = lineXEnd-lineXStart+1;
					
					
//					// The bounding rectangle
//					String lineColor = (r.strand==StrandedFeature.NEGATIVE) ? this.revStrandLineColor : this.fwStrandLineColor;
//					Element rectangle = doc.createElementNS(svgNS, "rect");
//					rectangle.setAttributeNS(null, "x", Integer.toString(lineXStart));
//					rectangle.setAttributeNS(null, "y", Integer.toString(lineYStart));
//					rectangle.setAttributeNS(null, "width", Integer.toString(lineWidth));
//					rectangle.setAttributeNS(null, "height", Integer.toString(lineHeight));
////					rectangle.setAttributeNS(null, "fill", "none");
//					rectangle.setAttributeNS(null, "fill", lineColor);
//					rectangle.setAttributeNS(null, "stroke", lineColor);
////					rectangle.setAttributeNS(null, "stroke-width", "1");
//					svgRoot.appendChild(rectangle);
					
					String lineColor = (r.strand==StrandedFeature.NEGATIVE) ? this.revStrandLineColor : this.fwStrandLineColor;
					Element line = doc.createElementNS(svgNS, "line");
					line.setAttributeNS(null, "x1", Double.toString(lineXStart));
					line.setAttributeNS(null, "y1", Integer.toString(rowCenter));
					line.setAttributeNS(null, "x2", Double.toString(lineXEnd));
					line.setAttributeNS(null, "y2", Integer.toString(rowCenter));
					line.setAttributeNS(null, "stroke", lineColor);
					line.setAttributeNS(null, "stroke-width", "1");
					svgRoot.appendChild(line);
					
					// Now go through individual CpGs
					for (int type = 0; type<=1; type++)
					{
						TreeSet<Integer> positions = (type==0) ? r.unmethPositions : r.methPositions;
						String fillColor = (type==0) ? this.unmethColor : this.methColor;
					
						for (Integer i : positions)
						{
							double pos = i.intValue();
							if (cpgPositionMap != null)
							{
								Double obj;
								obj = cpgPositionMap.get(i);
								if (obj == null) System.err.printf("Can't find read end position %.2f in cpgPosition map\n",pos);
								pos = obj.intValue();
							}

							double boxXStart = pos - firstCoord - (double)boxWidthHalf; 
							double boxXEnd = pos - firstCoord + (double)boxWidthHalf; 

							Element box = doc.createElementNS(svgNS, "rect");
							box.setAttributeNS(null, "x", Double.toString(boxXStart));
							box.setAttributeNS(null, "y", Integer.toString(boxYStart));
							box.setAttributeNS(null, "width", Integer.toString(this.boxWidth));
							box.setAttributeNS(null, "height", Integer.toString(boxHeight));
							box.setAttributeNS(null, "fill", fillColor);
							svgRoot.appendChild(box);
						}
					}
					

					

				

				}
				System.err.println();
			}
		}
		
//		
//		// Create the rectangle.
//		Element rectangle = doc.createElementNS(svgNS, "rect");
//		rectangle.setAttributeNS(null, "x", "10");
//		rectangle.setAttributeNS(null, "y", "20");
//		rectangle.setAttributeNS(null, "width", "100");
//		rectangle.setAttributeNS(null, "height", "50");
//		rectangle.setAttributeNS(null, "fill", "red");
//
//		// Attach the rectangle to the root 'svg' element.
//		svgRoot.appendChild(rectangle);
		
		
		return out;
	}
	
	

	private void addRead(MethylRead read) {
	
		boolean added = false;
		for (int l = 0; (l<MAX_LEVELS)&(!added); l++)
		{
			added = readFits(l, read);
			if (added) levels[l].add(read);
		}
		
		// if we got to the end, we don't have enough levels
		if (added == false)
		{
			System.err.printf("We need more than %d levels\n",MAX_LEVELS);
			(new Exception()).printStackTrace();
			System.exit(1);
		}
		
	}

	
	private boolean readFits(int level, MethylRead read) {

		boolean out = false;
		if (this.levels[level] == null)
		{
			// Nothing on this level yet
			out = true;
			this.levels[level] = new TreeSet<MethylRead>();
		}
		else
		{
			MethylRead last = this.levels[level].last();
			if (last == null)
			{
				System.err.printf("Why do we have a level %d with no elements??\n",level);
				(new Exception()).printStackTrace();
				System.exit(1);
			}
			out = read.startPos() > last.endPos();
		}
		
		return out;
	}

	public void writeDoc(Document doc, PrintWriter pw)
	{

        TransformerFactory tfactory = TransformerFactory.newInstance();
        Transformer serializer;
        try {
            serializer = tfactory.newTransformer();
            //Setup indenting to "pretty print"
            serializer.setOutputProperty(OutputKeys.INDENT, "yes");
            serializer.setOutputProperty("{http://xml.apache.org/xslt}indent-amount", "2");
            
            serializer.transform(new DOMSource(doc), new StreamResult(pw));
        } catch (TransformerException e) {
            // this is fatal, just dump the stack and throw a runtime exception
            e.printStackTrace();
            
            throw new RuntimeException(e);
        }
        
	}

	
	
}
