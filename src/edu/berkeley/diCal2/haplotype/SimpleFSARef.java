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
import gnu.trove.list.TIntList;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.map.TIntIntMap;
import gnu.trove.procedure.TIntIntProcedure;

import java.util.Arrays;

public class SimpleFSARef implements FSAReference {
	private final TIntIntMap localIdxMap;
	private final TIntList segList;
	
	private final byte[] alleleConfig;

	public SimpleFSARef (int[] alleleConfig, TIntIntMap localIdxMap) {
		this (SimpleFSARef.intArrayToByteArray(alleleConfig), localIdxMap);
	}
	
	public SimpleFSARef (byte[] alleleConfig, TIntIntMap localIdxMap) {
		this.alleleConfig = alleleConfig;
		this.localIdxMap = localIdxMap;

		// WARINING: for this to work, the range of localIdxMap has to be consecutive beginning at zero
		assert (hasConsecutiveRange  (localIdxMap, 0)) : "Range of index map is not consecutive. Maybe duplicate entries.";
		final int[] segArray = new int[localIdxMap.size()];
		localIdxMap.forEachEntry(new TIntIntProcedure() {
			public boolean execute(int arg0, int arg1) {
				segArray[arg1] = arg0;
				return true;
			}
		});
		
		this.segList = new TIntArrayList(segArray);
	}
	
	private static boolean hasConsecutiveRange (TIntIntMap map, int firstValue) {
		int[] range = map.values();
		Arrays.sort (range);
		for (int i=firstValue; i<range.length; i++) {
			if (i != range[i]) return false;
		}
		return true;
	}

	public boolean locusSegregating(int locus) {
		return localIdxMap.containsKey(locus);
	}

	public int getNonSegAllele(int locus) {
		return alleleConfig[locus];
	}

	public int getLocalIdx(int locus) {
		return localIdxMap.get(locus);
	}
	
	public int getSegLocus(int localIdx) {
		return segList.get(localIdx);
	}

	public int getNumLoci() {
		return alleleConfig.length;
	}

	public int getNumSegLoci() {
		return localIdxMap.size();
	}
	
	public int[] getSegSites () {
		return this.localIdxMap.keys();
	}

	public byte[] getAlleleConfig() {
		return this.alleleConfig;
	}
	
	public static int[] byteArrayToIntArray (byte[] byteArray) {
		int[] intArray = new int[byteArray.length];
		for (int i=0; i<byteArray.length; i++) {
			intArray[i] = byteArray[i];
		}
		return intArray;
	}
	
	public static byte[] intArrayToByteArray (int[] intArray) {
		byte[] byteArray = new byte[intArray.length];
		for (int i=0; i<intArray.length; i++) {
			byteArray[i] = (byte) intArray[i];
		}
		return byteArray;
	}
	
}
