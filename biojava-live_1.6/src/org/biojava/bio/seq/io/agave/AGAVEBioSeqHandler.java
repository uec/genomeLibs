/**
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
package org.biojava.bio.seq.io.agave;
import org.biojava.bio.BioException;
import org.biojava.bio.seq.DNATools;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.seq.StrandedFeature;
import org.biojava.bio.seq.impl.SimpleSequence;
import org.biojava.bio.symbol.Location;
import org.biojava.bio.symbol.SymbolList;
import org.biojava.utils.ChangeVetoException;
import org.xml.sax.Attributes;
import org.xml.sax.SAXException;

/**
 *
 * Handles the AGAVE &lt;bio_sequence&gt; element
 *
 * @author David Huen
 * @author Hanning Ni     Doubletwist Inc
 * @author Greg Cox
 */
public class AGAVEBioSeqHandler
               extends StAXFeatureHandler
               implements  AGAVEBioSeqCallbackItf, SequenceHandler
{
  public static final StAXHandlerFactory AGAVE_BIO_SEQ_HANDLER_FACTORY
    = new StAXHandlerFactory() {
    public StAXContentHandler getHandler(StAXFeatureHandler staxenv) {
      return new AGAVEBioSeqHandler(staxenv);
    }
  };

  //dna sequence,
  private SymbolList dna ;
  protected Sequence sequence ;

  AGAVEBioSeqHandler(StAXFeatureHandler staxenv) {
    // setup up environment stuff
    super( staxenv );
    featureListener = staxenv.featureListener;
    setHandlerCharacteristics("bio_sequence", true);

    // setup handlers
       // <db_id>
       super.addHandler(new ElementRecognizer.ByLocalName("db_id"),
         AGAVEDbIdPropHandler.AGAVE_DBID_PROP_HANDLER_FACTORY);
       // <note>
       super.addHandler(new ElementRecognizer.ByLocalName("note"),
         AGAVENotePropHandler.AGAVE_NOTE_PROP_HANDLER_FACTORY);
       // <gene>
       super.addHandler(new ElementRecognizer.ByLocalName("description"),
         AGAVEDescPropHandler.AGAVE_DESC_PROP_HANDLER_FACTORY);
       // <keyword>
       super.addHandler(new ElementRecognizer.ByLocalName("keyword"),
         AGAVEKeywordPropHandler.AGAVE_KEYWORD_PROP_HANDLER_FACTORY);
       // <sequence>
       super.addHandler(new ElementRecognizer.ByLocalName("sequence"),
         AGAVESeqPropHandler.AGAVE_SEQ_PROP_HANDLER_FACTORY);
       // <alt_ids>
       super.addHandler(new ElementRecognizer.ByLocalName("alt_ids"),
         AGAVEAltIdsPropHandler.AGAVE_ALT_IDS_PROP_HANDLER_FACTORY);
       // <xrefs>
       super.addHandler(new ElementRecognizer.ByLocalName("xrefs"),
         AGAVEXrefsPropHandler.AGAVE_XREFS_PROP_HANDLER_FACTORY);
       //<sequence_map>
       super.addHandler(new ElementRecognizer.ByLocalName("sequence_map"),
         AGAVESeqMapHandler.AGAVE_SEQ_MAP_HANDLER_FACTORY);
       //<map_location>
       super.addHandler(new ElementRecognizer.ByLocalName("map_location"),
         AGAVEMapLocationPropHandler.AGAVE_MAP_LOCATION_PROP_HANDLER_FACTORY);
       //<classification>
       super.addHandler(new ElementRecognizer.ByLocalName("classification"),
         AGAVEClassificationHandler.AGAVE_CLASSIFICATION_HANDLER_FACTORY);
  }

  public void reportStrand(StrandedFeature.Strand strand)
  {
    // obtains strand from elements that are in the know.
    ((StrandedFeature.Template) featureTemplate).strand = strand;
  }
  public void reportFeature(Location loc)
  {
    // obtains strand from elements that are in the know.
    ((StrandedFeature.Template) featureTemplate).location = loc;
  }


  public void reportDna(String dna_seq)
  {
    try{
       StringBuffer sb = new StringBuffer()  ;
       for( int i = 0 ; i < dna_seq.length(); i++)
       {
           char c = dna_seq.charAt(i) ;
           if( c != ' '  && c != '\n' && c!= '\t')
              sb.append( c  );
       }
       dna = DNATools.createDNA( sb.substring(0) );
     }catch(Exception e){ e.printStackTrace() ; }
  }

  public void startElementHandler(
                String nsURI,
                String localName,
                String qName,
                Attributes attrs)
                throws SAXException
  {
      try{
      featureListener.startSequence();
      boolean forFeature = false ;
      setProperty( "element_id",  attrs.getValue("element_id") , forFeature) ;
      setProperty( "sequence_id",  attrs.getValue("sequence_id") , forFeature) ;
      setProperty( "seq_length",  attrs.getValue("seq_length") , forFeature) ;
      setProperty( "molecule_type",  attrs.getValue("molecule_type") , forFeature) ;
      setProperty( "organism_name",  attrs.getValue("organism_name"), forFeature ) ;
      setProperty( "taxon_id",  attrs.getValue("taxon_id") ,   forFeature) ;
      setProperty( "clone_id",  attrs.getValue("clone_id"), forFeature ) ;
      setProperty( "clone_library",  attrs.getValue("clone_library"), forFeature ) ;
      setProperty( "chromosome",  attrs.getValue("chromosome") , forFeature) ;
      setProperty( "map_position",  attrs.getValue("map_position") , forFeature) ;
      setProperty( "ec_number",  attrs.getValue("ec_number"), forFeature ) ;
      setProperty( "create_date",  attrs.getValue("create_date") ,  forFeature) ;
      setProperty( "update_date",  attrs.getValue("update_date") , forFeature) ;
      }catch(Exception e)
      {
          throw new SAXException( e.getMessage() ) ;
      }
  }


   /*
   protected Feature.Template createTemplate() {
    // create Gene Template for this
    StrandedFeature.Template st = new StrandedFeature.Template();

    // assume feature set to describe a transcript
    st.type = "bio_sequence";
    st.strand = StrandedFeature.UNKNOWN;
    // set up annotation bundle
    st.annotation = new SmallAnnotation();
    st.location = new  Location.EmptyLocation();


    if( staxenv != null )
        staxenv. subFeatures .add( this ) ;

    return st;
  }*/


  public void endElementHandler(
                String nsURI,
                String localName,
                String qName,
                StAXContentHandler handler)
                throws SAXException
  {
       //super.endElement( nsURI, localName, qName, handler);
      //create sequence
      try{
       if( dna == null )
          throw new SAXException("dna sequence need offered for creating sequence object");

       sequence  = new SimpleSequence( dna, " ", "simple_sequence " , annot ) ;
       if( featureTemplate == null )
          throw new SAXException("feature template is null ") ;

       //Feature feature = sequence.createFeature(  featureTemplate ) ;
       //realizeSubFeatures( feature ) ;
       addFeatureToSequence(sequence) ;
       appendToTop(sequence, staxenv) ;
       featureListener.endSequence();
       }catch(BioException be){
            throw new SAXException( "bio exception" ) ;
       }catch(ChangeVetoException cve){
            throw new SAXException("change veto exception") ;
       }catch(Exception e)
       {
          throw new SAXException( e.getMessage() ) ;
       }
  }
   private void appendToTop(Sequence sequence, StAXFeatureHandler staxenv)
    {
       if (staxenv instanceof AGAVEContigCallbackItf)
       {
             ((AGAVEContigCallbackItf) staxenv).reportSequence( sequence );
             return;
       }
       if (staxenv instanceof AGAVECallbackItf)
       {
             ((AGAVECallbackItf) staxenv).reportSequence( sequence );
              return;
       }
       appendToTop(sequence, staxenv.staxenv );
    }

}

