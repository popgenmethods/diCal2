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

public class ConstDemography extends Demography {
	
	public ConstDemography (int numDemes, double migrationRate) {
		this (getUniformSizes(numDemes), migrationRate);
	}
	
	public ConstDemography (double[] popSizes, double migrationRate) {
		this (popSizes, RateMatrixTools.getSymmetricMigrationMatrix(migrationRate, popSizes.length));
	}
	
	public ConstDemography (double[] popSizes, double[][] migrationMatrix) {
		
		assert (popSizes.length == migrationMatrix.length);
		// just one long interval
		this.epochList = new RealPartition.Interval[] {new RealPartition.Interval (0, Double.POSITIVE_INFINITY)};

		// just everybody by themselves
		this.treePartitionList = new ArrayList<ArrayList<TreeSet<Integer>>>();
		// partition in long epoch
		ArrayList<TreeSet<Integer>> thePartition = new ArrayList<TreeSet<Integer>>();
		for (int i=0; i<popSizes.length; i++) {
			thePartition.add (new TreeSet<Integer>(Arrays.asList(new Integer[] {new Integer (i)})));
		}
		this.treePartitionList.add (thePartition);

		// just the one migration matrix forever
		this.migrationMatrixList = new ArrayList<double[][]>();
		this.migrationMatrixList.add (migrationMatrix);  
		
		// same pop sizes for all
		this.popSizes = new ArrayList<double[]>();
		this.popSizes.add (popSizes);
		
		// and no pulse migration
		this.pulseMigrationMatrixList = new ArrayList<double[][]>();
		this.pulseMigrationMatrixList.add (null);
	}

	private static double[] getUniformSizes (int numDemes) {
		double[] popSizes = new double[numDemes];
		for (int g=0; g<popSizes.length; g++) popSizes[g] = 1d/numDemes;
		return popSizes;
	}
	
}
