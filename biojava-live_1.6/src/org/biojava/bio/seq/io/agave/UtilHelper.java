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
package org.biojava.bio.seq.io.agave;
import org.biojava.bio.Annotation;

/**
 * Just make the property follow the common case
 * @author Hanning Ni   Doubletwist Inc
 */
public class UtilHelper {
   /**
    * inhibit the getProperty(key) of Annotation from throw exception when
    *  key does not exist.
    */
    public static Object getProperty(Annotation annot, String key)
    {
        if( annot == null || key == null) return null ;
        try{
           Object ob = annot.getProperty(key) ;
            return ob ;
        }catch( java.util.NoSuchElementException e)
        {
            return null ;
        }
    }
}