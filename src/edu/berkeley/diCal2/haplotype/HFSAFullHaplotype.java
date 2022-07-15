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


public class HFSAFullHaplotype extends HFSAHaplotype implements FSAFullHaplotype {
	public HFSAFullHaplotype(int[] alleleConfig) {
		super(alleleConfig);
	}

	public int getSpecifiedAllele(int locus) {
		return getAllele(locus);
	}
	
	public static class HFSAFullHapFactory implements FSAHapFactory<HFSAFullHaplotype> {
		public HFSAFullHaplotype getHaplotype(int[] alleleConfig) {
			return new HFSAFullHaplotype(alleleConfig.clone());
		}
	}
	
	public static class HFSAFullHapConfig extends GeneralFSAConfig<HFSAFullHaplotype> {
		public HFSAFullHapConfig(int numLoci, int numAlleles) {
			super(numLoci, numAlleles);
		}
		
		public void updateCachedItems(HFSAFullHaplotype haplotype, int adjustment) { }
		
		public static class HFSAFullHapConfigConstructor implements FSAConfigConstructor<HFSAFullHaplotype, HFSAFullHapConfig> {
			public HFSAFullHapConfig newFSAConfig(int numLoci, int numAlleles) {
				return new HFSAFullHapConfig(numLoci, numAlleles);
			}
		}
	}
}
