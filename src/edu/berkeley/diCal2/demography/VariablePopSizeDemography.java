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

package edu.berkeley.diCal2.demography;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.TreeSet;

import edu.berkeley.diCal2.utility.RateMatrixTools;
import edu.berkeley.diCal2.utility.RealPartition;

public class VariablePopSizeDemography extends Demography {
	public VariablePopSizeDemography (double[] popSizes, double[] changeTimes) {
		
		assert (popSizes.length-1 == changeTimes.length);
		
		// create intervals
		this.epochList = new RealPartition.Interval[changeTimes.length+1];
		double lastTime = 0d;
		for (int i=0; i<this.epochList.length; i++) {
			if (i < changeTimes.length) {
				this.epochList[i] = new RealPartition.Interval (lastTime, changeTimes[i]);
				// change last time
				lastTime = changeTimes[i];
			}
			else {
				this.epochList[i] = new RealPartition.Interval (lastTime, Double.POSITIVE_INFINITY);
			}
		}

		// just one partition for all times
		
		this.treePartitionList = new ArrayList<ArrayList<TreeSet<Integer>>>();
		// partition in long epoch
		for (int i=0; i<this.epochList.length; i++) {
			// make a singleton with one
			ArrayList<TreeSet<Integer>> thePartition = new ArrayList<TreeSet<Integer>>();
			thePartition.add (new TreeSet<Integer>(Arrays.asList(new Integer[] {new Integer (0)})));
			this.treePartitionList.add (thePartition);
		}

		// all three times the same migration matrix
		this.migrationMatrixList = new ArrayList<double[][]>();
		for (int i=0; i<this.epochList.length; i++) {
			// just take the symmetric guy
			this.migrationMatrixList.add (RateMatrixTools.getSymmetricMigrationMatrix (0d, 1));  
		}
		
		// same pop sizes for all
		this.popSizes = new ArrayList<double[]>();
		for (int i=0; i<this.epochList.length; i++) {
			// make an array with one entry
			double[] tmp = new double[1];
			tmp[0] = popSizes[i];
			this.popSizes.add (tmp);
		}
		
		// and some pulse migration
		this.pulseMigrationMatrixList = new ArrayList<double[][]>();
		for (int i=0; i<this.epochList.length; i++) {
			// no pulse migrations ar all
			this.pulseMigrationMatrixList.add(null);
		}
	}

}
