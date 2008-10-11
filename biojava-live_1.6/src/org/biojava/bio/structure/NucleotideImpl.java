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
 * Created on 05.03.2004
 * @author Andreas Prlic
 *
 */
package org.biojava.bio.structure;
import org.biojava.bio.structure.io.PDBParseException;


/**
 *
 *  A nucleotide group is almost the same as a Hetatm group. 
 *  @see HetatomImpl
 *  @see AminoAcidImpl
 * @author Andreas Prlic
 * @since 1.4
 * @version %I% %G%
 */
public class NucleotideImpl 
    extends HetatomImpl 
    implements Group			   
{

    /** this is a "nucleotide", a special occurance of a Hetatom. */
    public static final String type = "nucleotide";
   
    /*
     * inherits most from Hetero and has just a few extensions.
     */
    public NucleotideImpl() {
	super();

    }

    public String getType(){ return type;}

    	
    public String toString(){
		
	String str = "PDB: "+ pdb_name + " " + pdb_code +  " "+ pdb_flag;
	if (pdb_flag) {
	    str = str + "atoms: "+atoms.size();
	}
	return str ;
		
    }


    public Object clone(){
	NucleotideImpl n = new NucleotideImpl();

	n.setPDBFlag(has3D());
	n.setPDBCode(getPDBCode());
	try {
	    n.setPDBName(getPDBName());
	} catch (PDBParseException e) {
	    e.printStackTrace();
	}

	
	// copy the atoms
	for (int i=0;i<atoms.size();i++){
	    Atom atom = (Atom)atoms.get(i);
	    n.addAtom((Atom)atom.clone());
	}
	return n;
    }
}
