 /*
    This file is part of diCal2.

    diCal2 is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    diCal2 is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with diCal2.  If not, see <http://www.gnu.org/licenses/>.
  */

package edu.berkeley.diCal2.haplotype;

import edu.berkeley.diCal2.haplotype.FSAXTypeConfig.FSAReference;
import edu.berkeley.diCal2.utility.CollectionFormat;

public class HFSAXHaplotype implements FSAXHaplotype {
	protected final byte[] localAlleleConfig;
	protected final FSAReference fsaReference;

	private final int hashCodeCache;
	
	public HFSAXHaplotype(FSAReference fsaReference, byte[] localAlleleConfig) {
		this.localAlleleConfig = localAlleleConfig;
		this.fsaReference = fsaReference;
		
		this.hashCodeCache = getHashCode(fsaReference, localAlleleConfig);
	}

	public FSAReference getReference() {
		return fsaReference;
	}
	
	public int getAllele(int locus) {
		return fsaReference.locusSegregating(locus) ? localAlleleConfig[fsaReference.getLocalIdx(locus)] : fsaReference.getNonSegAllele(locus);
	}
	
	public int getAlleleByIdx(int locusIdx) {
		return localAlleleConfig[locusIdx];
	}
	
	public int getNumLoci() {
		return fsaReference.getNumLoci();
	}
	
	public int[] getRawXAlleleConfig() {
		return SimpleFSARef.byteArrayToIntArray(this.localAlleleConfig);
	}
	
	public int[] getRawAlleleConfig() {
		int[] rawAlleleConfig = new int[getNumLoci()];
		for (int l = 0; l < getNumLoci(); l++) {
			rawAlleleConfig[l] = getAllele(l);
		}
		
		return rawAlleleConfig;
	}
	
	public String toString() {
		return CollectionFormat.formatArray(SimpleFSARef.byteArrayToIntArray(this.localAlleleConfig), " ", "(", ")");
	}
	
	public int hashCode() {
		return hashCodeCache;
	}
	
	public static int getHashCode(FSAReference fsaReference, byte[] alleleConfig) {
		int result = 17;
		
		for (int l = 0; l < alleleConfig.length; l++) {
			result  = 37 * result + (alleleConfig[l] + 1);
		}

		result  = 37 * result + (System.identityHashCode(fsaReference) + 1);
		
		return result;
	}
	
	public final boolean equals(Object o) {
		
		if (o == this) return true;
		if (hashCodeCache != o.hashCode()) return false;
		if (!(o instanceof HFSAXHaplotype)) return false;
		
		HFSAXHaplotype toCompare = (HFSAXHaplotype) o;
		
		if (fsaReference != toCompare.fsaReference) return false;
		for (int l = 0; l < localAlleleConfig.length; l++) {
			if (localAlleleConfig[l] != toCompare.localAlleleConfig[l]) return false;
		}
		
		return true;
	}
	
	public static class HFSAXHapFactory implements FSAXHapFactory<HFSAXHaplotype> {
		public HFSAXHaplotype getHaplotype(FSAReference fsaReference, int[] alleleConfig) {
			return new HFSAXHaplotype(fsaReference, SimpleFSARef.intArrayToByteArray(alleleConfig));
		}
	}
	
	public static class HFSAXHapConfig extends GeneralFSAXConfig<HFSAXHaplotype> {
		public HFSAXHapConfig(int numLoci, int numAlleles, FSAReference fsaReference) {
			super(numLoci, numAlleles, fsaReference);
		}
		
		public void updateCachedItems(HFSAXHaplotype haplotype, int adjustment) { }
		
		public static class HFSAXHapConfigConstructor implements FSAXConfigConstructor<HFSAXHaplotype, HFSAXHapConfig> {
			public HFSAXHapConfig newFSAXConfig(int numLoci, int numAlleles, FSAReference fsaReference) {
				return new HFSAXHapConfig(numLoci, numAlleles, fsaReference);
			}
		}
	}
}
