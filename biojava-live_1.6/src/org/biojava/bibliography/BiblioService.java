// BiblioService.java
//
//    senger@ebi.ac.uk
//    March 2001
//

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
package org.biojava.bibliography;

/**
 * <p>
 * It represents a service dealing with the bibliographic
 * resources. It may identify a program name or some other
 * (usually) electronic agent who provides the cited resource.
 * </p>
 *
 * @author <A HREF="mailto:senger@ebi.ac.uk">Martin Senger</A>
 * @version $Id: BiblioService.java 2812 2003-07-16 16:01:11Z mrp $
 * @since 1.3
 */

public class BiblioService
    extends BiblioProvider {

  /**
   * The name of a bibliographic service.
   */
    public String name;
}
