/*
 *                  BioJava development code
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
 * Created on Dec 28, 2005
 *
 */
package org.biojava.utils.io;

import java.io.EOFException;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.util.Enumeration;
import java.util.jar.JarEntry;
import java.util.jar.JarFile;
import java.util.zip.GZIPInputStream;
import java.util.zip.ZipEntry;
import java.util.zip.ZipFile;


/** A class that provides an InputStream from a File. The file can be compressed or uncompressed. 
 *  
 * Currently supported
 * compressions:
 * <ul>
 * <li>Gzip (extension .gz)</li>
 * <li>Zip (extension .zip) in this case a stream to the first entry in the zip file is returned </li> 
 * <li>Jar (extension .jar) same as .Zip; only stream to first entry is returned </li>
 * <li>Z (extension .Z) compressed using the unix compress command </li>
 * <li>for any other extension, no compression is assumed </li>
 * </ul>
 * 
 * 
 * @author Andreas Prlic
 * @since 1.5
 * @version %I% %G%
 *
 */
public class InputStreamProvider {
    
    /**
     * The magic number found at the start of a GZIP stream.
     */
    public static final int GZIP_MAGIC = 0x1f8b;
    
    
    public InputStreamProvider() {
        super();
        
    }
    
    /** get an InputStream for this file 
     * 
     * @param pathToFile the path of the file.
     * @return an InputStream for the file located at the path.
     * @throws IOException
     */
    public InputStream getInputStream(String pathToFile)
    throws IOException
    {
        File f = new File(pathToFile);
        return getInputStream(f);
    }
        
    
    /** open the file and read the magic number from the beginning
     * this is used to determine the compression type
     * 
     * @param in an input stream to read from
     * @return the magic number
     * @throws IOException
     */
    private int getMagicNumber(InputStream in) 
    throws IOException {
              
        
        int t = in.read();
        if (t < 0) throw new EOFException("Failed to read magic number");
        int magic = (t & 0xff) << 8;
        t = in.read();
        if (t < 0) throw new EOFException("Failed to read magic number");
        magic += t & 0xff;
        
        return magic;
    }
    
    
    
    /** get an InputStream for the file
     * 
     * @param f a File
     * @return an InputStream for the file
     * @throws IOException
     */
    public InputStream getInputStream(File f) 
    throws IOException
    {
        // use the magic numbers to determine the compression type, 
        // use file extension only as 2nd choice 
        
        int magic = 0;
        
        try {
            InputStream test = new FileInputStream(f);
            magic = getMagicNumber(test);
            test.close();
        } catch (Exception e){
            e.printStackTrace();
        };
        
        InputStream inputStream = null;
       
        String fileName = f.getName();
        
        if (magic == UncompressInputStream.LZW_MAGIC ) {
            // a Z compressed file
            return openCompressedFile(f);
        }
        
        else if (magic == GZIP_MAGIC ) {
            return openGZIPFile(f); 
        }
        
        else if ( fileName.endsWith(".gz")) {
            return openGZIPFile(f);
        } 
        
        else if ( fileName.endsWith(".zip")){
          
            ZipFile zipfile = new ZipFile(f);
            
            // stream to first entry is returned ...
            ZipEntry entry;
            Enumeration e = zipfile.entries();
            if ( e.hasMoreElements()){
                entry = (ZipEntry) e.nextElement();
                inputStream = zipfile.getInputStream(entry);
                
            } else {
                throw new IOException ("Zip file has no entries");
            }
            
        } 
        
        else if ( fileName.endsWith(".jar")) {

            JarFile jarFile = new JarFile(f);
           
            // stream to first entry is returned
            JarEntry entry;
            Enumeration e = jarFile.entries();
            if ( e.hasMoreElements()){
                entry = (JarEntry) e.nextElement();
               
                inputStream = jarFile.getInputStream(entry);
            } else {
                throw new IOException ("Jar file has no entries");
            }
        } 
        
        else if ( fileName.endsWith(".Z")) {
            // unix compressed 
            return openCompressedFile(f);
            
        }
        
        else {
           
            // no particular extension found, assume that it is an uncompressed file
            inputStream = new FileInputStream(f);
        }
        
        return inputStream;
    }
    

    
    private InputStream openCompressedFile(File f)
    throws IOException{
        
        InputStream is           = new FileInputStream(f);
        InputStream inputStream =  new UncompressInputStream(is);
        return inputStream;
    }
    
    private InputStream openGZIPFile(File f) 
    throws IOException{
        
        FileInputStream is      = new FileInputStream(f);
        InputStream inputStream = new GZIPInputStream(is);
        return inputStream;
    }
   
}
