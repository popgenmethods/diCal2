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


public class HFSAXFullHaplotype extends HFSAXHaplotype implements FSAXFullHaplotype {
	public HFSAXFullHaplotype(FSAReference fsaReference, int[] localAlleleConfig) {
		this (fsaReference, SimpleFSARef.intArrayToByteArray(localAlleleConfig));
	}

	public HFSAXFullHaplotype(FSAReference fsaReference, byte[] localAlleleConfig) {
		super(fsaReference, localAlleleConfig);
	}

	public int getSpecifiedAllele(int locus) {
		return getAllele(locus);
	}
	
	public static class HFSAXFullHapFactory implements FSAXHapFactory<HFSAXFullHaplotype> {
		public HFSAXFullHaplotype getHaplotype(FSAReference fsaReference, int[] alleleConfig) {
			return new HFSAXFullHaplotype(fsaReference, alleleConfig);
		}
	}
	
	public static class HFSAXFullHapConfig extends GeneralFSAXConfig<HFSAXFullHaplotype> {
		public HFSAXFullHapConfig(int numLoci, int numAlleles, FSAReference fsaReference) {
			super(numLoci, numAlleles, fsaReference);
		}
		
		public void updateCachedItems(HFSAXFullHaplotype haplotype, int adjustment) { }
		
		public static class HFSAXFullHapConfigConstructor implements FSAXConfigConstructor<HFSAXFullHaplotype, HFSAXFullHapConfig> {
			public HFSAXFullHapConfig newFSAXConfig(int numLoci, int numAlleles, FSAReference fsaReference) {
				return new HFSAXFullHapConfig(numLoci, numAlleles, fsaReference);
			}
		}
	}
	
	/// wrapper calls to be able to emulate an empty haplotype
	public static class HFSAXFullHaplotypeShell implements FSAXFullHaplotype {

		// this is the haplotype
		private HFSAXFullHaplotype hap;
		
		// we need a unique number
		private static int numHaps = 0;
		private final int hapID;
		private int getAndIncreaseNumHaps() {
			return ++numHaps;
		}
		public int getHapID () {
			return this.hapID;
		}
		
		public HFSAXFullHaplotypeShell (HFSAXFullHaplotype hap) {
			// store the given haplotype
			this.hap = hap;
			// give it a unique number
			this.hapID = getAndIncreaseNumHaps();
		}
		
		// if we don't want the allele data anymore
		public void clear () {
			this.hap = null;
		}
		
		@Override
		public int getAlleleByIdx(int locusIdx) {
			assert (this.hap != null);
			return this.hap.getAlleleByIdx (locusIdx);
		}

		@Override
		public int[] getRawXAlleleConfig() {
			assert (this.hap != null);
			return this.hap.getRawXAlleleConfig();
		}

		@Override
		public FSAReference getReference() {
			assert (this.hap != null);
			return this.hap.getReference();
		}

		@Override
		public int getNumLoci() {
			assert (this.hap != null);
			return this.hap.getNumLoci();
		}

		@Override
		public int getAllele(int locus) {
			assert (this.hap != null);
			return this.hap.getAllele(locus);
		}

		@Override
		public int[] getRawAlleleConfig() {
			assert (this.hap != null);
			return this.hap.getRawAlleleConfig();
		}

		@Override
		public int getSpecifiedAllele(int locus) {
			assert (this.hap != null);
			return this.hap.getSpecifiedAllele(locus);
		}
		
		public String toString () {
			if (this.hap != null) {
				return this.hap.toString();
			}
			else {
				return "[hap deleted " + this.getHapID() + " ]";
			}
		}
	}
}
