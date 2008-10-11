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
 * Created on 16.03.2004
 * @author Andreas Prlic
 *
 *
 */
package org.biojava.bio.structure.io;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.net.URL;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.zip.GZIPOutputStream;

import org.biojava.bio.structure.AminoAcid;
import org.biojava.bio.structure.Chain;
import org.biojava.bio.structure.Group;
import org.biojava.bio.structure.GroupIterator;
import org.biojava.bio.structure.Structure;
import org.biojava.utils.io.InputStreamProvider;


/**
 * <p>
 *  The wrapper class for parsing a PDB file.
 *  </p>
 *  
 *  
 *  <p>
 *  Several flags can be set for this class
 *  <ul>
 *  <li> {@link #setParseCAOnly} - parse only the Atom records for C-alpha atoms (default:false)</li>
 * <li> {@link #setParseSecStruc} - a flag if the secondary structure information from the PDB file (author's assignment) should be parsed.
 *      If true the assignment can be accessed through {@link AminoAcid}.getSecStruc(); (default:false)</li>
 * <li> {@link #setAlignSeqRes(boolean)} - should the AminoAcid sequences from the SEQRES
 *      and ATOM records of a PDB file be aligned? (default:true)</li>    
 * <li> {@link #setAutoFetch(boolean)} - if the PDB file can not be found locally, should it be fetched
 *  from the EBI - ftp server? (default:false)</li>
 *  </ul>
 *  </p>
 * 
 *
 *
 *<h2>Example</h2>
 * <p>
 * Q: How can I get a Structure object from a PDB file?
 * </p>
 * <p>
 * A:
 * <pre>
 public {@link Structure} loadStructure(String pathToPDBFile){
		{@link PDBFileReader} pdbreader = new {@link PDBFileReader}();

		{@link Structure} structure = null;
		try{
			structure = pdbreader.getStructure(pathToPDBFile);
			System.out.println(structure);
		} catch (IOException e) {
			e.printStackTrace();
		}
		return structure;
	}
 </pre>
 *
 * Access PDB files from a directory, take care of compressed PDB files
 * <pre>
 * public {@link Structure} loadStructureById() {
		String path = "/path/to/PDB/directory/";

		{@link PDBFileReader} pdbreader = new {@link PDBFileReader}();
		pdbreader.setPath(path);
		{@link Structure} structure = null;
		try {
			structure = pdbreader.getStructureById("5pti");
		} catch (IOException e){
			e.printStackTrace();
		}
		return structure;
		
	}
	</pre>
 * 
 *
 * @author Andreas Prlic
 *
 */
public class PDBFileReader implements StructureIOFile {

// a list of big pdb files for testing
//  "1htq",
//  "1c2w",
//  "1ffk",
//  "1giy",
//  "1j5a",
//  "1jj2",
//  "1jzx",
//  "1jzy",
//  "1jzz",
//  "1k01",
//  "1k73",
//  "1k8a",
//  "1k9m",
//  "1kc8",
//  "1kd1",
//  "1kqs",
//  "1m1k",
//  "1m90",
//  "1mkz",
//  "1ml5",
//  "1n8r",
  
  
  
    
	String path                     ;
	List<String> extensions            ;
	boolean parseSecStruc;
	boolean autoFetch;
	boolean parseCAOnly;
    boolean alignSeqRes;
  
    
	public static void main(String[] args){
		String filename =  "/path/to/PDBFile.pdb" ;

		// also see the demos
		
		PDBFileReader pdbreader = new PDBFileReader();
		pdbreader.setParseSecStruc(true);
		pdbreader.setAlignSeqRes(true);
		pdbreader.setParseCAOnly(false);
		pdbreader.setAutoFetch(true);

		try{
			Structure struc = pdbreader.getStructure(filename);
			System.out.println(struc);

			GroupIterator gi = new GroupIterator(struc);
			while (gi.hasNext()){
				Group g = (Group) gi.next();
				Chain  c = g.getParent();
				if ( g instanceof AminoAcid ){
					AminoAcid aa = (AminoAcid)g;
					Map<String,String> sec = aa.getSecStruc();
					System.out.println(c.getName() + " " + g + " " + sec);
				}

			}

		} catch (Exception e) {
			e.printStackTrace();
		}
	}


	public PDBFileReader() {
		extensions    = new ArrayList<String>();
		path = "" ;
		extensions.add(".ent");
		extensions.add(".pdb");
		extensions.add(".ent.gz");
		extensions.add(".pdb.gz");
		extensions.add(".ent.Z");
		extensions.add(".pdb.Z");
		parseSecStruc = false;
		autoFetch     = false;
        parseCAOnly   = false;
        alignSeqRes   = true;
	}


    /** return the flag if only the CA atoms should be parsed
     * 
     * @return flag if CA only should be read
     */
	public boolean isParseCAOnly() {
        return parseCAOnly;
    }

	/** only the CA atoms should be parsed from the PDB file
     * 
     * @param parseCAOnly
	 */
    public void setParseCAOnly(boolean parseCAOnly) {
        this.parseCAOnly = parseCAOnly;
    }

    /** get the flag if the SEQRES and ATOM amino acids are going to be aligned
     * 
     * @return flag
     */
    public boolean isAlignSeqRes() {
        return alignSeqRes;
    }


    /** set the flag if the SEQRES and ATOM amino acids should be aligned and linked
     * 
     * @param alignSeqRes
     */
    public void setAlignSeqRes(boolean alignSeqRes) {
        this.alignSeqRes = alignSeqRes;
    }


    /** should the parser to fetch missing PDB files from the EBI FTP server automatically?
	 *  default is false
	 * @return flag
	 */
	public boolean isAutoFetch() {
		return autoFetch;
	}

	/** tell the parser to fetch missing PDB files from the EBI FTP server automatically.
	 * 
	 * default is false. If true, new PDB files will be automatically stored in the Path and gzip compressed.
	 * 
	 * @param autoFetch
	 */
	public void setAutoFetch(boolean autoFetch) {
		this.autoFetch = autoFetch;
	}

	/* A flag to tell the parser to parse the Author's secondary structure assignment from the file
	 * default is set to false, i.e. do NOT parse.
	 * @param parseSecStruc
	 */
	public boolean isParseSecStruc() {
		return parseSecStruc;
	}


	/*  A flag to tell the parser to parse the Author's secondary structure assignment from the file
	 *  
	 */
	 
	public void setParseSecStruc(boolean parseSecStruc) {
		this.parseSecStruc = parseSecStruc;
	}


	/** directory where to find PDB files */
	public void setPath(String p){
		path = p ;
	}

	/**
	 * Returns the path value.
	 * @return a String representing the path value
	 * @see #setPath
	 *
	 */
	public String getPath() {
		return path ;
	}

	/** define supported file extensions 
	 * compressed extensions .Z,.gz do not need to be specified
	 * they are dealt with automatically.

	 */
	public void addExtension(String s){
		//System.out.println("add Extension "+s);
		extensions.add(s);
	}




	/** try to find the file in the filesystem and return a filestream in order to parse it 
	 * rules how to find file
	 * - first check: if file is in path specified by PDBpath
	 * - secnd check: if not found check in PDBpath/xy/ where xy is second and third char of PDBcode.
	 * if autoFetch is set it will try to download missing PDB files automatically.
	 */

	private InputStream getInputStream(String pdbId) 
	throws IOException
	{
		//System.out.println("checking file");

		// compression formats supported
		// this has been moved to InputStreamProvider ...

		//String[] str = {".gz",".zip",".Z"};
		//ArrayList  compressions = new ArrayList( Arrays.asList( str ) ); 

		InputStream inputStream =null;

		String pdbFile = null ;
		File f = null ;

		// this are the possible PDB file names...
		String fpath = path+"/"+pdbId;
		String ppath = path +"/pdb"+pdbId;

		String[] paths = new String[]{fpath,ppath};

		for ( int p=0;p<paths.length;p++ ){
			String testpath = paths[p];
			//System.out.println(testpath);
			for (int i=0 ; i<extensions.size();i++){
				String ex = (String)extensions.get(i) ;
				//System.out.println("PDBFileReader testing: "+testpath+ex);
				f = new File(testpath+ex) ;

				if ( f.exists()) {
					//System.out.println("found!");
					pdbFile = testpath+ex ;

					InputStreamProvider isp = new InputStreamProvider();

					inputStream = isp.getInputStream(pdbFile);
					break;
				}

				if ( pdbFile != null) break;        
			}
		}

		if ( pdbFile == null ) {
			if ( autoFetch) 
				return downloadAndGetInputStream(pdbId);

			String message = "no structure with PDB code " + pdbId + " found!" ;
			throw new IOException (message);
		}

		return inputStream ;
	}


	private File downloadPDB(String pdbId){
	
		File tempFile = new File(path+"/"+pdbId+".pdb.gz");		
		File pdbHome = new File(path);
		
		if ( ! pdbHome.canWrite() ){
			System.err.println("can not write to " + pdbHome);
			return null;
		}
		
		String ftp = String.format("ftp://ftp.ebi.ac.uk/pub/databases/msd/pdb_uncompressed/pdb%s.ent", pdbId.toLowerCase()); 

		System.out.println("Fetching " + ftp); 
		try {
			URL url = new URL(ftp);
			InputStream conn = url.openStream();			

			// prepare destination
			System.out.println("writing to " + tempFile);
			
		
			FileOutputStream outPut = new FileOutputStream(tempFile);
			GZIPOutputStream gzOutPut = new GZIPOutputStream(outPut);
			PrintWriter pw = new PrintWriter(gzOutPut);
			
			BufferedReader fileBuffer = new BufferedReader(new InputStreamReader(conn));
			String line;
			while ((line = fileBuffer.readLine()) != null) {				
				pw.println(line);
			}
			pw.flush();
			pw.close();
			outPut.close();
			conn.close();
		} catch (Exception e){
			e.printStackTrace();
			return null; 
		}
		return tempFile;
	}

	private InputStream downloadAndGetInputStream(String pdbId) 
	throws IOException{
		//PDBURLReader reader = new PDBURLReader();
		//Structure s = reader.getStructureById(pdbId);
		File tmp = downloadPDB(pdbId);
		if ( tmp != null ) {
			InputStreamProvider prov = new InputStreamProvider();
			return prov.getInputStream(tmp);
			

		} else {
			throw new IOException("could not find PDB " + pdbId + " in file system and also could not download");
		}


	}




	/** load a structure from local file system and return a PDBStructure object 

	 * @param pdbId  a String specifying the id value (PDB code)
	 * @return the Structure object
	 * @throws IOException ...
	 */
	public Structure getStructureById(String pdbId) 
	throws IOException
	{

	    
		InputStream inStream = getInputStream(pdbId);

		PDBFileParser pdbpars = new PDBFileParser();
		pdbpars.setParseSecStruc(parseSecStruc);
        pdbpars.setAlignSeqRes(alignSeqRes);
        
        pdbpars.setParseCAOnly(parseCAOnly);
        
        
		Structure struc = pdbpars.parsePDBFile(inStream) ;
		return struc ;
	}

	/** opens filename, parses it and returns
	 * aStructure object .
	 * @param filename  a String
	 * @return the Structure object
	 * @throws IOException ...
	 */
	public Structure getStructure(String filename) 
	throws IOException
	{
		File f = new File(filename);
		return getStructure(f);

	}

	/** opens filename, parses it and returns a Structure object
	 * 
	 * @param filename a File object
	 * @return the Structure object
	 * @throws IOException ...
	 */
	public Structure getStructure(File filename) throws IOException {

		InputStreamProvider isp = new InputStreamProvider();

		InputStream inStream = isp.getInputStream(filename);

		PDBFileParser pdbpars = new PDBFileParser();
		pdbpars.setParseSecStruc(parseSecStruc);
        pdbpars.setAlignSeqRes(alignSeqRes);
        pdbpars.setParseCAOnly(parseCAOnly);
		Structure struc = pdbpars.parsePDBFile(inStream) ;
		return struc ;

	}

}
