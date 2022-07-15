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

package edu.berkeley.diCal2.utility;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;
import java.io.Reader;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.Random;

import edu.berkeley.diCal2.csd.auxiliary.DelayedRandom;


public class VisitPermuations {
	public static class Permutation {
		private final int[] permutation;
		
		public Permutation(int[] permutation) {
			this.permutation = permutation;
		}
		
		public <T> void permute(ArrayList<T> srcPermute, ArrayList<T> dstPermute) {
			for (int pIdx = 0; pIdx < permutation.length; pIdx++) {
				dstPermute.set(pIdx, srcPermute.get(permutation[pIdx]));
			}
		}
		
		public <T> ArrayList<T> permute(ArrayList<T> srcPermute) {
			ArrayList<T> dstPermute = new ArrayList<T>(permutation.length);
			for (int pIdx = 0; pIdx < permutation.length; pIdx++) {
				dstPermute.add(srcPermute.get(permutation[pIdx]));
			}
			
			return dstPermute;
		}
		
		public static List<Permutation> readPermutationList(Reader reader) throws NumberFormatException, IOException {
			List<Permutation> permutationList = new ArrayList<Permutation>();

			BufferedReader bufferedReader = new BufferedReader(reader);

			String nextLine = null;
			while ((nextLine = bufferedReader.readLine()) != null) {
				String[] preArray = nextLine.split(" ");
				// also split at tabs
				List<String> permutationArray = new ArrayList<String>();
				for (String curr : preArray) {
					for (String piece : curr.split("\t")) {
						permutationArray.add(piece);
					}
				}
				
				int[] rawPermutation = new int[permutationArray.size()];
				
				for (int pIdx = 0; pIdx < permutationArray.size(); pIdx++) {
					rawPermutation[pIdx] = Integer.parseInt(permutationArray.get(pIdx));
				}
				
				permutationList.add(new Permutation(rawPermutation));
			}
			
			return permutationList;
		}
		
		public static void writePermutationList(List<Permutation> permutationList, OutputStream oStream) {
			PrintStream printStream = new PrintStream(oStream);
			
			for (Permutation permutation : permutationList) {
				StringBuilder permString = new StringBuilder();
				
				for (int p : permutation.permutation) {
					permString.append(p + " ");
				}

				printStream.println(permString);
			}
		}
		
		public static List<Permutation> generatePermutationList(int numPermutations, int n, DelayedRandom rGen) {
			return generatePermutationList(numPermutations, n, n, rGen);
		}
		
		public static List<Permutation> generatePermutationList(int numPermutations, int n, int k, DelayedRandom rGen) {
			List<Permutation> permutationList = new ArrayList<Permutation>(numPermutations);
			
			ArrayList<Integer> shuffleList = new ArrayList<Integer>(n);
			for (int i = 0; i < n; i++) shuffleList.add(i);
			
			for (int p = 0; p < numPermutations; p++) {
				Collections.shuffle(shuffleList, rGen.getInternalRandom());
				
				int[] permutation = new int[k];
				for (int i = 0; i < k; i++) permutation[i] = shuffleList.get(i);
				
				permutationList.add(new Permutation(permutation));
			}
			
			return permutationList;
		}
	}
	
	public interface PermutationVisitor<T> {
		public boolean visitPermutation(List<T> permutation);
	}
	
	public static <T> void visitPermutations(int numPermuations, Collection<T> toPermute, PermutationVisitor<T> permutationVisitor, Random rGen) {
		ArrayList<T> shuffleList = new ArrayList<T>(toPermute);
		
		for (int p = 0; p < numPermuations; p++) {
			Collections.shuffle(shuffleList, rGen);
			
			permutationVisitor.visitPermutation(shuffleList);
		}
	}
	
	public static <T> void visitPermutations(Collection<T> toPermute, PermutationVisitor<T> permutationVisitor, List<Permutation> permutationList) {
		ArrayList<T> srcList = new ArrayList<T>(toPermute);
		ArrayList<T> dstList = new ArrayList<T>(toPermute);
		
		for (Permutation permutation : permutationList) {
			permutation.permute(srcList, dstList);
			
			permutationVisitor.visitPermutation(dstList);
		}
	}
}
