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

import java.io.PrintStream;

public class RateMatrixTools {

	public static boolean isSquare (double[][] matrix) {
		int numRows = matrix.length;
		for (int i=0; i<numRows; i++) {
			if (matrix[i].length != numRows) return false;
		}
		// we made it through, so everything fine
		return true;
	}
	
	/// check whether matrix is square, non-negative off-diagonal entries, and the sum of the off diagonal is the negative of the diagonal
	public static boolean isProperRateMatrix (double[][] matrix, double EPSILON) {
		double sum;
		int i,j;
		// first everything is fine
		boolean result = true;
		
		// get the column dimension
		int colDim = matrix.length;
		
		// go through and see whether everything is fine
		for (i=0; i<colDim; i++) {
			result &= (colDim == matrix[i].length);
			sum = 0.0;
			for (j=0; j<colDim; j++) {
				if (i != j) {
					// positive rate?
					result &= (matrix[i][j] >= 0.0);
					sum += matrix[i][j];
				}
			}
			// is the negative off-diagonal sum equal to diagonal entry? (proper rate matrix?)
			result &= (Math.abs(sum + matrix[i][i]) < EPSILON);
		}
		
		// return it?
		return result;
	}

	/// check whether a given matrix 
	public static boolean isProperStochaticMatrix (double[][] matrix, double EPSILON) {
		// just use the rate matrix thing for matrix - Id
		// copy it and subtract the identity
		double[][] localMatrix = new double[matrix.length][matrix.length];
		for (int i=0; i<matrix.length;  i++) {
			assert (matrix.length == matrix[i].length);
			for (int j=0; j<matrix[i].length; j++) {
				if ((i != j) && (matrix[i][j] > 1d + EPSILON)) return false;
				localMatrix[i][j] = matrix[i][j];
				// subtract one, if we need to
				if (i == j)  localMatrix[i][i] -= 1;
			}
		}
		// now check the matrix
		return isProperRateMatrix(localMatrix, EPSILON);
	}
	
	/// check whether all entries are zero
	public static boolean allEntriesZero(double[][] matrix, double EPSILON) {
		boolean result = true;
		for (double[] values : matrix) {
			for (double value : values) {
				result &= (Math.abs (value) < EPSILON);
			}
		}
		return result;
	}
	
	/// dump matrix
	public static void dump (double[][] matrix, PrintStream outStream) {
		for (double[] values : matrix) {
			for (double value : values) {
				outStream.print (value + "\t");
			}
			outStream.println();
		}
	}
	
	// this shall help thou
	public static double[][] getSymmetricMigrationMatrix (double betweenDemeMigRate, int numDemes) {

		assert (numDemes > 0);

		double[][] demeMigRates = new double[numDemes][numDemes];

		if (numDemes == 1) {
			// nothing flows nowhere
			demeMigRates = new double[][] {{0d}};
		}
		else {
			// here is some flowing
			for (int srcDemeIdx = 0; srcDemeIdx < numDemes; srcDemeIdx++) {
				for (int dstDemeIdx = 0; dstDemeIdx < numDemes; dstDemeIdx++) {
					demeMigRates[srcDemeIdx][dstDemeIdx] = (srcDemeIdx == dstDemeIdx) ? - (numDemes - 1) * betweenDemeMigRate : betweenDemeMigRate;
				}
			}
		}

		return demeMigRates;
	}


}
