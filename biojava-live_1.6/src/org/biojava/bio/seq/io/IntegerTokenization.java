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


package org.biojava.bio.seq.io;

import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.bio.symbol.IntegerAlphabet;
import org.biojava.bio.symbol.Symbol;

/**
 * @author Thomas Down
 */
public class IntegerTokenization extends WordTokenization {
    public IntegerTokenization() {
	super(IntegerAlphabet.getInstance());
    }

    public Symbol parseToken(String seq)
        throws IllegalSymbolException
    {
	try {
	    int i = Integer.parseInt(seq);
	    return ((IntegerAlphabet) getAlphabet()).getSymbol(i);
	} catch (NumberFormatException ex) {
	    throw new IllegalSymbolException(ex, "Couldn't parse " + seq);
	}
    }

    public String tokenizeSymbol(Symbol sym)
        throws IllegalSymbolException
    {
	getAlphabet().validate(sym);
	return sym.getName();
    }
}
