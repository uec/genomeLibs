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

package org.biojava.utils.xml;

import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

/**
 * Implementation of XMLWriter which emits nicely formatted documents
 * to a PrintWriter.
 *
 * @author Thomas Down
 * @since 1.3
 */

public class PrettyXMLWriter implements XMLWriter {
    private int indentUnit = 2;

    private PrintWriter writer;
    private boolean isOpeningTag = false;
    private boolean afterNewline = true;
    private int indent = 0;
    
    private Map namespacePrefixes = new HashMap();
    private int namespaceSeed = 0;
    private LinkedList namespaceBindings = new LinkedList();
    private List namespacesDeclared = new ArrayList();

    public PrettyXMLWriter(PrintWriter writer) {
        this.writer = writer;
    }

    public void declareNamespace(String nsURI, String prefixHint) 
        throws IOException
    {
        if (!namespacePrefixes.containsKey(nsURI)) {
            if (isOpeningTag == true) {
                String prefix = allocPrefix(nsURI);
                attribute("xmlns:" + prefix, nsURI);
            } else {
                namespacesDeclared.add(nsURI);
            }
        }
    }
    
    private void handleDeclaredNamespaces() 
        throws IOException
    {
        if (namespacesDeclared.size() == 0) {
            for (Iterator nsi = namespacesDeclared.iterator(); nsi.hasNext(); ) {
                String nsURI = (String) nsi.next();
                if (!namespacePrefixes.containsKey(nsURI)) {
                    String prefix = allocPrefix(nsURI);
                    attribute("xmlns:" + prefix, nsURI);
                }
            }
            namespacesDeclared.clear();
        }
    }
    
    protected void writeIndent()
        throws IOException
    {
        for (int i = 0; i < indent * indentUnit; ++i) {
            writer.write(' ');
        }
    }

    private void _openTag()
        throws IOException
    {
        if (isOpeningTag) {
            writer.println('>');
            afterNewline = true;
        }
        if (afterNewline) {
            writeIndent();
        }
        indent++;
        isOpeningTag = true;
        afterNewline = false;
        namespaceBindings.add(null);
    }
    
    private String allocPrefix(String nsURI) {
        String prefix = "ns" + (++namespaceSeed);
        namespacePrefixes.put(nsURI, prefix);
        List bindings = (List) namespaceBindings.getLast();
        if (bindings == null) {
            bindings = new ArrayList();
            namespaceBindings.removeLast();
            namespaceBindings.add(bindings);
        }
        bindings.add(nsURI);
        return prefix;
    }
    
    public void openTag(String nsURI, String localName)
        throws IOException
    {
    	if (nsURI == null || nsURI.length() == 0)
    	{
    		throw new IOException("Invalid namespace URI: "+nsURI);
    	}
        _openTag();
        boolean alloced = false;
        String prefix = (String) namespacePrefixes.get(nsURI);
        if (prefix == null) {
            prefix = allocPrefix(nsURI);
            alloced = true;
        }
        writer.print('<');
        writer.print(prefix);
        writer.print(':');
        writer.print(localName);
        if (alloced) {
            attribute("xmlns:" + prefix, nsURI);
        }
        handleDeclaredNamespaces();
    }
    
    public void openTag(String qName)
        throws IOException
    {
        _openTag();
        writer.print('<');
        writer.print(qName);
        handleDeclaredNamespaces();
    }

    public void attribute(String nsURI, String localName, String value)
        throws IOException
    {
        if (! isOpeningTag) {
            throw new IOException("attributes must follow an openTag");
        }

        String prefix = (String) namespacePrefixes.get(nsURI);
        if (prefix == null) {
            prefix = allocPrefix(nsURI);
            attribute("xmlns:" + prefix, nsURI);
        }
        
        writer.print(' ');
        writer.print(prefix);
        writer.print(':');
        writer.print(localName);
        writer.print("=\"");
        printAttributeValue(value);
        writer.print('"');
    }
    
    public void attribute(String qName, String value)
        throws IOException
    {
        if (! isOpeningTag) {
            throw new IOException("attributes must follow an openTag");
        }

        writer.print(' ');
        writer.print(qName);
        writer.print("=\"");
        printAttributeValue(value);
        writer.print('"');
    }

    private void _closeTag() {
        isOpeningTag = false;
        afterNewline = true;
        List hereBindings = (List) namespaceBindings.removeLast();
        if (hereBindings != null) {
            for (Iterator bi = hereBindings.iterator(); bi.hasNext(); ) {
                namespacePrefixes.remove(bi.next());
            }
        }
    }
    
    public void closeTag(String nsURI, String localName)
        throws IOException
    {
        String prefix = (String) namespacePrefixes.get(nsURI);
        if (prefix == null) {
            throw new IOException("Assertion failed: unknown namespace when closing tag");
        }
        indent--;

        if (isOpeningTag) {
            writer.println(" />");
        } else {
            if (afterNewline) {
                writeIndent();
            }
            writer.print("</");
            writer.print(prefix);
            writer.print(':');
            writer.print(localName);
            writer.println('>');
        }
        _closeTag();
    }
    
    public void closeTag(String qName) 
        throws IOException
    {
        indent--;

        if (isOpeningTag) {
            writer.println(" />");
        } else {
            if (afterNewline) {
                writeIndent();
            }
            writer.print("</");
            writer.print(qName);
            writer.println('>');
        }
        _closeTag();
    }

    public void println(String data)
        throws IOException
    {
	if (isOpeningTag) {
	    writer.println('>');
	    isOpeningTag = false;
	}
	printChars(data);
	writer.println();
	afterNewline = true;
    }

    public void print(String data)
        throws IOException
    {
	if (isOpeningTag) {
	    writer.print('>');
	    isOpeningTag = false;
	}
	printChars(data);
	afterNewline = false;
    }

    public void printRaw(String data)
        throws IOException
    {
	writer.println(data);
    }

    protected void printChars(String data) 
        throws IOException
    {
	if (data == null) {
	    printChars("null");
	    return;
	}

	for (int pos = 0; pos < data.length(); ++pos) {
	    char c = data.charAt(pos);
	    if (c == '<' || c == '>' || c == '&') {
		numericalEntity(c);
	    } else {
		writer.write(c);
	    }
	}
    }

    protected void printAttributeValue(String data) 
        throws IOException
    {
	if (data == null) {
	    printAttributeValue("null");
	    return;
	}

	for (int pos = 0; pos < data.length(); ++pos) {
	    char c = data.charAt(pos);
	    if (c == '<' || c == '>' || c == '&' || c == '"') {
		numericalEntity(c);
	    } else {
		writer.write(c);
	    }
	}
    }

    protected void numericalEntity(char c)
        throws IOException
    {
	writer.print("&#");
	writer.print((int) c);
	writer.print(';');
    }
    
    public void close()
        throws IOException
    {
        writer.close();
    }
}
