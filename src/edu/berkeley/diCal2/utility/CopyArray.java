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

import java.util.Arrays;
import java.util.Collections;

public class CopyArray {
	public static double[] nCopies(int numCopies, double d) {
		double[] nCopyArray = new double[numCopies];
		Arrays.fill(nCopyArray, d);
		
		return nCopyArray;
	}
	
	public static int[] nCopies(int numCopies, int d) {
		int[] nCopyArray = new int[numCopies];
		Arrays.fill(nCopyArray, d);
		
		return nCopyArray;
	}

	public static boolean[] nCopies(int numCopies, boolean b) {
		boolean[] nCopyArray = new boolean[numCopies];
		Arrays.fill(nCopyArray, b);
		
		return nCopyArray;
	}
	
	public static <O> O[] nCopies(int numCopies, O obj, O[] objArr) {
		return Collections.nCopies(numCopies, obj).toArray(objArr);
	}

}
