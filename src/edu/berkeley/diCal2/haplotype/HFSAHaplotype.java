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

import edu.berkeley.diCal2.utility.CollectionFormat;

public class HFSAHaplotype implements FSAHaplotype {
	protected final int[] alleleConfig;
	private final int hashCodeCache;
	
	public HFSAHaplotype(int[] alleleConfig) {
		this.alleleConfig = alleleConfig;
		this.hashCodeCache = getHashCode(alleleConfig);
	}

	public int getAllele(int locus) {
		return alleleConfig[locus];
	}
	
	public int getNumLoci() {
		return alleleConfig.length;
	}
	
	public int[] getRawAlleleConfig() {
		return alleleConfig.clone();
	}
	
	public String toString() {
		return CollectionFormat.formatArray(alleleConfig, " ", "(", ")");
	}
	
	public int hashCode() {
		return hashCodeCache;
	}
	
	public static int getHashCode(int[] alleleConfig) {
		int result = 17;
		
		for (int l = 0; l < alleleConfig.length; l++) {
			result  = 37 * result + (alleleConfig[l] + 1);
		}
		
		return result;
	}
	
	public final boolean equals(Object o) {
		
		if (o == this) return true;
		if (hashCodeCache != o.hashCode()) return false;
		if (!(o instanceof HFSAHaplotype)) return false;
		
		HFSAHaplotype toCompare = (HFSAHaplotype) o;
		for (int l = 0; l < alleleConfig.length; l++) {
			if (alleleConfig[l] != toCompare.alleleConfig[l]) return false;
		}
		
		return true;
	}
	
	public static class HFSAHapFactory implements FSAHapFactory<HFSAHaplotype> {
		public HFSAHaplotype getHaplotype(int[] alleleConfig) {
			return new HFSAHaplotype(alleleConfig.clone());
		}
	}
	
	public static class HFSAHapConfig extends GeneralFSAConfig<HFSAHaplotype> {
		public HFSAHapConfig(int numLoci, int numAlleles) {
			super(numLoci, numAlleles);
		}
		
		public void updateCachedItems(HFSAHaplotype haplotype, int adjustment) { }
		
		public static class HFSAHapConfigConstructor implements FSAConfigConstructor<HFSAHaplotype, HFSAHapConfig> {
			public HFSAHapConfig newFSAConfig(int numLoci, int numAlleles) {
				return new HFSAHapConfig(numLoci, numAlleles);
			}
		}
	}
}
