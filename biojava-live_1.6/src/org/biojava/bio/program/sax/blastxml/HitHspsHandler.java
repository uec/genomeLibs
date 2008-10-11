/*
 *  BioJava development code This code may be freely distributed and modified
 *  under the terms of the GNU Lesser General Public Licence. This should be
 *  distributed with the code. If you do not have a copy, see:
 *  http://www.gnu.org/copyleft/lesser.html Copyright for this code is held
 *  jointly by the individual authors. These should be listed in
 *
 *@author    doc comments. For more information on the BioJava project and its
 *      aims, or to join the biojava-l mailing list, visit the home page at:
 *      http://www.biojava.org/
 */
package org.biojava.bio.program.sax.blastxml;

import org.biojava.bio.seq.io.game.ElementRecognizer;
import org.biojava.utils.stax.StAXContentHandler;
import org.xml.sax.Attributes;
import org.xml.sax.SAXException;
import org.xml.sax.helpers.AttributesImpl;

/**
 * @author David Huen
 */
class HitHspsHandler
    extends StAXFeatureHandler
{
    // create static factory class that makes an instance
    // of this class.
    public final static StAXHandlerFactory HIT_HSPS_HANDLER_FACTORY
             =
        new StAXHandlerFactory() {
            public StAXContentHandler getHandler(StAXFeatureHandler staxenv) {
                return new HitHspsHandler(staxenv);
            }
        };

    // constructor
    public HitHspsHandler(StAXFeatureHandler staxenv)
    {
        super(staxenv);
//        System.out.println("HitHspsHandler staxenv " + staxenv);
        // delegate handling of <Hsp> to its own class
        super.addHandler(new ElementRecognizer.ByLocalName("Hsp"),
            HspHandler.HSP_HANDLER_FACTORY);
    }

    public void startElementHandler(
            String nsURI,
            String localName,
            String qName,
            Attributes attrs)
             throws SAXException
    {
        // generate start of <biojava:HSPCollection>
        staxenv.listener.startElement(biojavaUri, "HSPCollection", biojavaUri + ":HSPCollection", new AttributesImpl());        
    }

    public void endElementHandler(
            String nsURI,
            String localName,
            String qName,
            StAXContentHandler handler)
             throws SAXException
    {
        // generate end of <biojava:HSPCollection>
        staxenv.listener.endElement(biojavaUri, "HSPCollection", biojavaUri + ":HSPCollection");
    }
}
