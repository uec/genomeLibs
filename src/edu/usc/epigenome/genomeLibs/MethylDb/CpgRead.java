/**
 * 
 */
package edu.usc.epigenome.genomeLibs.MethylDb;

/**
 * @author benb
 *
 */
public class CpgRead {


	// Object vars
	public int readId = 0;
	
	public short cRead = 0;
	public short cReadNonconversionFilt = 0;
	public short tRead = 0;
	public short agRead = 0;
	
	public short nextBaseGread = 0;
	protected char nextBaseUpperCase = '0';

	
	/**
	 * @param readId
	 */
	public CpgRead(int readId) {
		super();
		this.readId = readId;
		
	}

	/**
	 * @param readId
	 * @param cRead
	 * @param cReadNonconversionFilt
	 * @param tRead
	 * @param agRead
	 */
	public CpgRead(int readId, short cRead, short cReadNonconversionFilt,
			short tRead, short agRead, short nextBaseGread, char nextBaseUpperCase) {
		super();
		this.readId = readId;
		this.cRead = cRead;
		this.cReadNonconversionFilt = cReadNonconversionFilt;
		this.tRead = tRead;
		this.agRead = agRead;
		this.nextBaseGread = nextBaseGread;
		this.nextBaseUpperCase = nextBaseUpperCase;
	}


	
}
