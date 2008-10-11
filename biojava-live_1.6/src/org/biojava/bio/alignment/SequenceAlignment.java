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

/*
 * Created on 2005-08-03
 *
 */
package org.biojava.bio.alignment;

import java.util.List;
import org.biojava.bio.seq.SequenceIterator;
import org.biojava.bio.seq.db.SequenceDB;
import org.biojava.bio.symbol.Alignment;
import org.biojava.bio.symbol.SymbolList;

/**
 * This Interface provides methods for the alignment of biosequences.
 * 
 * @author Andreas Dr&auml;ger
 * @author Mark Schreiber
 * 
 */
public abstract class SequenceAlignment {
  /**
   * @return a string representation of the alignment
   * @throws BioException
   */
  public abstract String getAlignmentString() throws Exception;

  /**
   * @param source
   *          a SequenceIterator containing a set of sequences to be aligned
   *          with
   * @param subjectDB
   *          the SequenceDB containing another set of sequences.
   * @return a list containing the results of all single alignments performed by
   *         this method.
   * @throws NoSuchElementException
   * @throws Exception
   */
  public abstract List alignAll(SequenceIterator source, SequenceDB subjectDB)
      throws Exception;

  /**
   * Performs a pairwise sequence alignment of the two given sequences.
   * 
   * @param query
   * @param subject
   * @return score of the alignment or the distance.
   * @throws Exception
   */
  public abstract double pairwiseAlignment(SymbolList query, SymbolList subject)
      throws Exception;

  /**
   * This method also performs a sequence alignment of the two given sequences
   * but it returns an Alignment object instead of the score.
   * 
   * @param query
   * @param subject
   * @return Alignment
   */
  public abstract Alignment getAlignment(SymbolList query, SymbolList subject)
      throws Exception;

  /**
   * This method provides a BLAST-like formated alignment from the given
   * <code>String</code>s, in which the sequence coordinates and the
   * information "Query" or "Target", respectively is added to each line. Each
   * line contains 60 sequence characters including the gap symbols plus the
   * meta information. There is one white line between two pairs of sequences.
   * 
   * @param queryName
   *          name of the query sequence
   * @param targetName
   *          name of the target sequence
   * @param align
   *          a <code>String</code>-array, where the index 0 is the query
   *          sequence and index 1 the target sequence (for instance
   *          <code>new String[] {myQuerySequence.seqString(), myTargetSequence.seqString()}</code>)
   * @param path
   *          the "path", that means a String containing white spaces and pipe
   *          ("|") symbols, which makes matches visible. The two strings in
   *          <code>align</code> have to have the same length and also the
   *          same length than this <code>path</code>.
   * @param queryStart
   *          the start position in the query, where the alignment starts. For
   *          example zero for normal Needleman-Wunsch-Alignments.
   * @param queryEnd
   *          the end position, that means the sequence coordinate, which is the
   *          last symbol of the query sequence. Counting starts at zero!
   * @param queryLength
   *          The length of the query sequence without gaps.
   * @param targetStart
   *          These are all the same for the target. Have a look at these above.
   * @param targetEnd
   * @param targetLength
   * @param editdistance
   * @param time
   *          The time in milliseconds, which was needed to generate the
   *          alignment.
   * @return formatierten String.
   */
  public static String formatOutput(String queryName, String targetName,
      String[] align, String path, int queryStart, int queryEnd,
      long queryLength, int targetStart, int targetEnd, long targetLength,
      double editdistance, long time) {
    String output = System.getProperty("line.separator") + " Time (ms):\t"
        + time + System.getProperty("line.separator") + " Length:\t"
        + align[0].length() + System.getProperty("line.separator");
    output += "  Score:\t" + (-1) * editdistance
        + System.getProperty("line.separator");
    output += "  Query:\t" + queryName + ",\tLength:\t" + queryLength
        + System.getProperty("line.separator");
    output += "  Target:\t" + targetName + ",\tLength:\t" + targetLength
        + System.getProperty("line.separator")
        + System.getProperty("line.separator");

    int currline = Math.min(60, align[0].length()), i, j, k, l;
    // counts the absolute position within the String
    String space = "  ", kspace = "", jspace = "";
    for (k = 0; k < new Integer(Math.max(queryEnd, targetEnd)).toString()
        .length(); k++)
      space += " ";
    for (k = new Integer(queryStart + 1).toString().length(); k <= new Integer(
        Math.max(queryEnd, targetEnd)).toString().length(); k++)
      kspace += " ";
    for (k = new Integer(targetStart + 1).toString().length(); k <= new Integer(
        Math.max(queryEnd, targetEnd)).toString().length(); k++)
      jspace += " ";

    i = k = queryStart;
    j = l = targetStart;
    output += System.getProperty("line.separator") + "Query:\t" + kspace
        + (k + 1) + " ";
    for (i = currline - Math.min(60, align[0].length()); i < currline; i++) {
      if ((align[0].charAt(i) != '-') && (align[0].charAt(i) != '~'))
        k++;
      if ((align[1].charAt(i) != '-') && (align[1].charAt(i) != '~'))
        j++;
    }
    output += align[0].substring(0, currline) + " " + k;
    output += " " + System.getProperty("line.separator") + "        " + space
        + path.substring(0, currline);
    output += " " + System.getProperty("line.separator") + "Target:\t" + jspace
        + (l + 1) + " " + align[1].substring(0, currline) + " " + j + " "
        + System.getProperty("line.separator");

    for (; currline + 60 < path.length(); currline += 60) {
      l = Math.min(j + 1, targetEnd);
      kspace = jspace = "";
      for (int n = new Integer(k + 1).toString().length() - 1; n < new Integer(
          Math.max(queryEnd, targetEnd)).toString().length(); n++)
        kspace += " ";
      for (int n = new Integer(j).toString().length() - 1; n < new Integer(Math
          .max(queryEnd, targetEnd)).toString().length(); n++)
        jspace += " ";
      output += " " + System.getProperty("line.separator") + "Query:\t"
          + kspace + Math.min(k + 1, queryEnd) + " ";
      for (i = currline; i < currline + 60; i++) {
        if ((align[0].charAt(i) != '-') && (align[0].charAt(i) != '~'))
          k++;
        if ((align[1].charAt(i) != '-') && (align[1].charAt(i) != '~'))
          j++;
      }
      output += align[0].substring(currline, currline + 60) + " " + k;
      output += " " + System.getProperty("line.separator") + "        " + space
          + path.substring(currline, currline + 60);
      output += " " + System.getProperty("line.separator") + "Target:\t"
          + jspace + l + " " + align[1].substring(currline, currline + 60)
          + " " + j + " " + System.getProperty("line.separator");
    }
    align[0] += " " + queryEnd;
    align[1] += " " + targetEnd;
    if (currline + 1 < path.length()) {
      kspace = jspace = "";
      for (int n = new Integer(k).toString().length() - 1; n < new Integer(Math
          .max(queryEnd, targetEnd)).toString().length(); n++)
        kspace += " ";
      for (int n = new Integer(j).toString().length() - 1; n < new Integer(Math
          .max(queryEnd, targetEnd)).toString().length(); n++)
        jspace += " ";
      output += " " + System.getProperty("line.separator") + "Query:\t"
          + kspace + Math.min(k + 1, queryEnd) + " "
          + align[0].substring(currline, align[0].length());
      output += " " + System.getProperty("line.separator") + "        " + space
          + path.substring(currline, path.length());
      output += " " + System.getProperty("line.separator") + "Target:\t"
          + jspace + Math.min(j + 1, targetEnd) + " "
          + align[1].substring(currline, align[1].length())
          + System.getProperty("line.separator");
    }

    return output += System.getProperty("line.separator");
  }

}
