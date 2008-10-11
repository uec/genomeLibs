/*
 *                    BioJava development code
 *
 * This code may be freely distributed and modified under the
 * terms of the GNU Lesser General Public Licence.  This should
 * be distributed with the code.  If you do not have a copy,
 * see:
 *
 *      http://www.gnu.org/copyleft/lesser.html
 *
 * Copyright for this code is held jointly by the individual
 * authors.  These should be listed in @author doc comments.
 *
 * For more information on the BioJava project and its aims,
 * or to join the biojava-l mailing list, visit the home page
 * at:
 *
 *      http://www.biojava.org/
 *
 */
package org.biojava.bio.program.sax;

import java.util.HashMap;
import java.util.StringTokenizer;

import org.xml.sax.SAXException;

/**
 * A Helper class for parsing summary lines of Blast-like
 * output. For example:
 *
 * Y00365 Chinese hamster high mobility group prote...   246  1e-62
 *
 * The number of tokens from the right-hand side to collect
 * is specified by the calling method. In the above case,
 * the 2 right-most tokens i.e. 246 and 1e-62 would be chosen.
 *
 * Primary author -
 *                 Simon Brocklehurst (CAT)
 * Other authors  -
 *                 Tim Dilks          (CAT)
 *                 Colin Hardman      (CAT)
 *                 Stuart Johnston    (CAT)
 *
 * Copyright 2000 Cambridge Antibody Technology Group plc.
 *
 *
 * This code released to the biojava project, May 2000
 * under the LGPL license.
 *
 * @author Cambridge Antibody Technology Group plc
 * @author Greg Cox
 * @version 0.1
 *
 */
final class NcbiBlastSummaryLineHelper implements SummaryLineHelperIF{

    public NcbiBlastSummaryLineHelper() {
    }

    /**
     * Takes a line such as
     *
     * L38477 Mus musculus (clone Clebp-1) high...   353  7e-95
     *
     * And parses it according to the following rules
     *
     * From the left, and tokenizing on white space, extracts the
     * first token (above, this would be "L38477") and places
     * it as a String in the object Buffer.
     *
     * From the right, and tokenizing on white space, looks for
     * a specified number of tokens, which it places as
     * Strings in the object map
     *
     * @param poLine	 -
     * @param poMap	 A HashMap of name-value pairs to be
     * be interpreted by the calling class. The first two
     * items in the map will be the HitId and the HitDescription.
     * Subsequent will be attribute name-values pairs such as
     * Score, E-value.
     */
    public void parse(String poLine, HashMap poMap,
		      BlastLikeVersionSupport poVersion) throws SAXException {


	int iGrab = 2;
	int iCount;
	StringBuffer oHitDescription = new StringBuffer();

	StringTokenizer oSt = new StringTokenizer(poLine);


	//ncbi-blast all flavours - two tokens:
	//first is score
	//next is Evalue


	iGrab = 2;


	//populate Map...
	iCount = oSt.countTokens() - iGrab - 1;
	//first token is the hit id
	poMap.put("hitId",oSt.nextToken());
	oHitDescription.setLength(0);

	for (int i = 0; i < iCount; i++) {
	    oHitDescription.append(oSt.nextToken());
	    oHitDescription.append(" ");
	}
	poMap.put("hitDescription",oHitDescription.substring(0));

	//now collect score and e-value
	poMap.put("score",oSt.nextToken());
	poMap.put("expectValue",oSt.nextToken());

    }

}

