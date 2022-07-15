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

public class ExampleDemography extends Demography {
	
	public static class ExampleDemography2 extends ExampleDemography {
		
		public ExampleDemography2() {
			super();
			
			this.expGrowthRates = new ArrayList<double[]>();
			for (int epoch = 0; epoch < this.epochList.length; epoch++) {
				this.expGrowthRates.add(new double[this.popSizes.get(epoch).length]);
			}
		}
		
	}
	
	// our example
	public ExampleDemography () {
		
		// first some partition of the real line to determine the epochs
		this.epochList = new RealPartition.Interval[] {
													new RealPartition.Interval (0,0.5),
													new RealPartition.Interval (0.5,1),
													new RealPartition.Interval (1,1.5),
													new RealPartition.Interval (1.5,Double.POSITIVE_INFINITY),
												};

		// now some tree structure to describe demes in each epoch
		this.treePartitionList = new ArrayList<ArrayList<TreeSet<Integer>>>();
		// first epoch
		ArrayList<TreeSet<Integer>> firstInterval = new ArrayList<TreeSet<Integer>>();
		firstInterval.add (new TreeSet<Integer>(Arrays.asList(new Integer[] {new Integer (0)})));
		firstInterval.add (new TreeSet<Integer>(Arrays.asList(new Integer[] {new Integer (1)})));
		firstInterval.add (new TreeSet<Integer>(Arrays.asList(new Integer[] {new Integer (2)})));
		firstInterval.add (new TreeSet<Integer>(Arrays.asList(new Integer[] {new Integer (3)})));
		this.treePartitionList.add (firstInterval);
		// second epoch
		ArrayList<TreeSet<Integer>> secondInterval = new ArrayList<TreeSet<Integer>>();
		secondInterval.add (new TreeSet<Integer>(Arrays.asList(new Integer[] {new Integer (0)})));
		secondInterval.add (new TreeSet<Integer>(Arrays.asList(new Integer[] {new Integer (1)})));
		secondInterval.add (new TreeSet<Integer>(Arrays.asList(new Integer[] {new Integer (2), new Integer (3)})));
		this.treePartitionList.add (secondInterval);
		// third epoch
		ArrayList<TreeSet<Integer>> thirdInterval = new ArrayList<TreeSet<Integer>>();
		thirdInterval.add (new TreeSet<Integer>(Arrays.asList(new Integer[] {new Integer (0)})));
		thirdInterval.add (new TreeSet<Integer>(Arrays.asList(new Integer[] {new Integer (1), new Integer (2), new Integer (3)})));
		this.treePartitionList.add (thirdInterval);
		// fourth
		ArrayList<TreeSet<Integer>> fourthInterval = new ArrayList<TreeSet<Integer>>();
		fourthInterval.add (new TreeSet<Integer>(Arrays.asList(new Integer[] {new Integer (0), new Integer (1), new Integer (2), new Integer (3)})));
		this.treePartitionList.add (fourthInterval);

		// now let's get some migration matrices and pop sizes
		// for now we just take symmetric migration among all demes that are present
		double migrationRate = 1d;
		migrationMatrixList = new ArrayList<double[][]>();
		// add all 4
		migrationMatrixList.add (RateMatrixTools.getSymmetricMigrationMatrix (migrationRate, treePartitionList.get(0).size()));  
		migrationMatrixList.add (RateMatrixTools.getSymmetricMigrationMatrix (migrationRate, treePartitionList.get(1).size()));  
		migrationMatrixList.add (RateMatrixTools.getSymmetricMigrationMatrix (migrationRate, treePartitionList.get(2).size()));  
		migrationMatrixList.add (RateMatrixTools.getSymmetricMigrationMatrix (migrationRate, treePartitionList.get(3).size()));  
		
		// the sizes of the demes
		popSizes = new ArrayList<double[]>();
		popSizes.add (new double[] {0.25, 0.25, 0.25, 0.5});
		popSizes.add (new double[] {1d/3, 1d/3, 1d/3});
		popSizes.add (new double[] {0.75, 0.5});
		popSizes.add (new double[] {1.0});
		
		// and no pulse migrations
		this.pulseMigrationMatrixList = new ArrayList<double[][]>(); 
		for (int i=0; i<this.epochList.length; i++) {
			this.pulseMigrationMatrixList.add (null);
		}

	}
	
}
