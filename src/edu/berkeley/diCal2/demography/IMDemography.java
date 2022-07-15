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

public class IMDemography extends Demography {
	public IMDemography (double splittingTime, double stopGeneFlowTime, double symmetricMigrationRate, double[] initialPopSizes, double[] geneFlowPopSizes, double ancestralPopSize) {
		
		assert (stopGeneFlowTime < splittingTime);
		// first some partition of the real line to determine the epochs
		this.epochList = new RealPartition.Interval[] {
													new RealPartition.Interval (0, stopGeneFlowTime),
													new RealPartition.Interval (stopGeneFlowTime, splittingTime),
													new RealPartition.Interval (splittingTime, Double.POSITIVE_INFINITY),
												};

		// now some tree structure to describe demes in each epoch
		this.treePartitionList = new ArrayList<ArrayList<TreeSet<Integer>>>();
		// first epoch
		ArrayList<TreeSet<Integer>> firstInterval = new ArrayList<TreeSet<Integer>>();
		firstInterval.add (new TreeSet<Integer>(Arrays.asList(new Integer[] {new Integer (0)})));
		firstInterval.add (new TreeSet<Integer>(Arrays.asList(new Integer[] {new Integer (1)})));
		this.treePartitionList.add (firstInterval);
		// second epoch
		ArrayList<TreeSet<Integer>> secondInterval = new ArrayList<TreeSet<Integer>>();
		secondInterval.add (new TreeSet<Integer>(Arrays.asList(new Integer[] {new Integer (0)})));
		secondInterval.add (new TreeSet<Integer>(Arrays.asList(new Integer[] {new Integer (1)})));
		this.treePartitionList.add (secondInterval);
		// third epoch
		ArrayList<TreeSet<Integer>> thirdInterval = new ArrayList<TreeSet<Integer>>();
		thirdInterval.add (new TreeSet<Integer>(Arrays.asList(new Integer[] {new Integer (0), new Integer (1)})));
		this.treePartitionList.add (thirdInterval);

		// now let's get some migration matrices and pop sizes
		// for now we just take symmetric migration among all demes that are present
		this.migrationMatrixList = new ArrayList<double[][]>();
		// add all 4
		this.migrationMatrixList.add (RateMatrixTools.getSymmetricMigrationMatrix (0.0, treePartitionList.get(0).size()));  
		this.migrationMatrixList.add (RateMatrixTools.getSymmetricMigrationMatrix (symmetricMigrationRate, treePartitionList.get(1).size()));  
		this.migrationMatrixList.add (RateMatrixTools.getSymmetricMigrationMatrix (symmetricMigrationRate, treePartitionList.get(2).size()));  
		
		assert(initialPopSizes.length == 2);
		assert(geneFlowPopSizes.length == 2);
		
		this.popSizes = new ArrayList<double[]>();
		this.popSizes.add (Arrays.copyOf(initialPopSizes, initialPopSizes.length));
		this.popSizes.add (Arrays.copyOf(geneFlowPopSizes, geneFlowPopSizes.length));
		this.popSizes.add (new double[] {ancestralPopSize});
		
		// no pulse migration yet
		this.pulseMigrationMatrixList = new ArrayList<double[][]>();
		for (int i=0; i<this.epochList.length; i++) {
			this.pulseMigrationMatrixList.add (null);
		}
	}
}
