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

package edu.berkeley.diCal2.csd;

import java.util.TreeMap;

public class SMCSDTransitionList {
	// [currIntervalDeme][prevIntervalDeme]
	private final TreeMap<Integer, double[][]> logRecoTransition = new TreeMap<Integer, double[][]>();
	private final TreeMap<Integer, double[]> logNoRecoTransition = new TreeMap<Integer, double[]>();

	public double[][] getLogReco (int index) {
		double[][] returnValue = logRecoTransition.get(index);
		assert (returnValue != null);
		return returnValue;
	}
	public void setLogReco (int index, double[][] value) {
		logRecoTransition.put (index, value);
	}

	public double[] getLogNoReco (int index) {
		double[] returnValue = logNoRecoTransition.get(index);
		assert (returnValue != null);
		return returnValue;
	}
	public void setLogNoReco (int index, double[] value) {
		logNoRecoTransition.put (index, value);
	}
}
