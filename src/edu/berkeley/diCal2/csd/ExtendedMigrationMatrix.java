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

import edu.berkeley.diCal2.utility.RateMatrixTools;
import Jama.Matrix;

public class ExtendedMigrationMatrix {

	/// construct an extended migration param set from old migration guy and a list of sample sizes
	public static double[][] computeExtendedMigrationMatrix (double[][] mMatrix, double[] absorptionRates) {
		// we know that migration things proper, so build around it
		
		// I guess it does not hurt to store the old one
		assert (RateMatrixTools.isProperRateMatrix (mMatrix, UberDemographyCore.ONELOCUS_EPSILON));
		
		int numDemes = mMatrix.length;
		assert (numDemes == absorptionRates.length);
		// make some range (for later use)
		int indexRange[] = new int[numDemes];
		for (int i=0; i<numDemes; i++) {
			indexRange[i] = i; 
		}

		// make some N matrix (Jama matrix somehow useful here)
		Matrix generalizedMigMatrix = new Matrix (2*numDemes, 2*numDemes, 0.0);
		// get the coalescent rates for migration, that is the migration rates divided by two
		Matrix migrationRates = new Matrix(mMatrix);
		
		generalizedMigMatrix.setMatrix (indexRange, indexRange, migrationRates);
		for (int i=0; i<numDemes; i++) {
			// subtract absorption rate from diagonal
			double tmp = generalizedMigMatrix.get (i, i) - absorptionRates[i];
			generalizedMigMatrix.set (i, i, tmp);
			// and add in the right place
			generalizedMigMatrix.set (i, i+numDemes, absorptionRates[i]);
		}
		
		// store the matrix
		double[][] extendedMigrationMatrix = generalizedMigMatrix.getArray();
		
		// just to be sure
		assert (RateMatrixTools.isProperRateMatrix (extendedMigrationMatrix, UberDemographyCore.ONELOCUS_EPSILON));
		
		// give it away now
		return extendedMigrationMatrix;
	}
}
