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
import org.biojava.bio.seq.Sequence;

/**
 * An interface that can be tested for by nested handlers
 * when trying to do a callback.
 * This one handles callbacks from nested elements that
 * determine strandedness of a nesting element.
 *
 *
 * @author Hanning Ni     Doubletwist Inc
 */
public interface AGAVEChromosomeCallbackItf
{
    public void reportSequence(Sequence sequence);
}

