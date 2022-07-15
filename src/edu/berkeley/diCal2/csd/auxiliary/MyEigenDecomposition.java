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

package edu.berkeley.diCal2.csd.auxiliary;

import org.netlib.lapack.DGEEV;
import org.netlib.util.intW;

import edu.berkeley.diCal2.utility.RateMatrixTools;
import Jama.Matrix;
import Jampack.Inv;
import Jampack.JampackException;
import Jampack.Plus;
import Jampack.Times;
import Jampack.Z;
import Jampack.Zmat;

public class MyEigenDecomposition {

	/// the precision
	private final double EPSILON;
	
	/// list of possibly complex eigenvectors
	final public Z[] lambda;
	
	/// let us also store the matrix of vectors
	final public Zmat V;
	
	/// and its inverse
	final public Zmat Vinv;
	
	/// list of the outer products of the eigenvectors and the row vectors of the inverse [order matches the eigenvalues]
	final public Zmat[] VW;

	final private double[][] origMatrix;
	
	// so we don't have to change the gold standard
	public MyEigenDecomposition(double[][] origMatrix) {
		this(origMatrix, 1e-10);
	}
	
	/// build an extended eigendecomposition object from a given matrix
	public MyEigenDecomposition(double[][] origMatrix, double EPSILON) {
		
		this.EPSILON = EPSILON;
		
		// just remember it
		this.origMatrix = origMatrix;
		
		// check quadraticness and get dimensions
		int matDim = this.origMatrix.length;
		for (double[] row : this.origMatrix) {
			assert (row.length == matDim);
		}

		// I guess we should make a copy of the original matrix, cause the eigenstuff might mess up the matrix
		double[][] matrix = new double[matDim][matDim];
		for (int i=0; i<matDim; i++) {
			for (int j=0; j<matDim; j++) {
				matrix[i][j] = this.origMatrix[i][j];
			}
		}
		
		// use jlapack for basic eigendecomposition

		// the sandboxes for the eigenstuff
		double[] realValues = new double[matDim];
		double[] imaginaryValues = new double[matDim];
		double[][] leftVectors = new double[matDim][matDim];
		double[][] rightVectors = new double[matDim][matDim];
		double[] sandBox = new double[5*matDim];
		intW returnInt = new intW(0);
		// calculate some eigenstuff (hopefully the correct call)
		DGEEV.DGEEV ("N", "V", matDim, matrix, realValues, imaginaryValues, leftVectors, rightVectors, sandBox, sandBox.length, returnInt);

		// make the diagonal guy D with eigenvalues
		lambda = new Z[matDim];
		for (int i=0; i<matDim; i++) {
			Z tmp = new Z (realValues[i], imaginaryValues[i]);
			lambda[i] = tmp;
		}

		
		// get the V matrix and the diagonal matrix of eigenvalues
		// pay special attention to potential complex values though
		// get the possibly complex eigenvectors into a suitable jampack matrix
		// we can actually initialize the matrix with the rightVectors as real part
		// and then just modify the complex vectors accordingly
		V = new Zmat (rightVectors);
		// go through index and adjust unreql vectors
		double realPart, imaginaryPart;
		// initialize the counter
		int idx = 0;
		while (idx<matDim) {
			
			// is the current value just real?
			// positive guy should be first, so we need just one-directional check
			if (lambda[idx].im > EPSILON) {
				// i and i+1 is unreal
				assert (Math.abs(lambda[idx].re - lambda[idx+1].re) < EPSILON);
				assert (Math.abs(lambda[idx].im + lambda[idx+1].im) < EPSILON);
				
				// modify the vectors at i and i+1
				for (int j=0; j<matDim; j++) {
					//get the values [they are stored in this way]
					realPart = V.get0(j,idx).re;
					assert (Math.abs (V.get0(j,idx).im) < EPSILON);
					imaginaryPart = V.get0(j,idx+1).re;
					assert (Math.abs (V.get0(j,idx+1).im) < EPSILON);
					
					// set the values right [hopefully we get it right]
					V.put0(j,idx, new Z (realPart, imaginaryPart));
					V.put0(j,idx+1, new Z (realPart, -imaginaryPart));
				}

				// increase i by two, since their was a pair of complex guys
				idx += 2;
			}
			else {
				// just one real guy, just increase index by one
				idx += 1;
			}
		}

		// check the eigenvalues and eigenvectors whether we can assume this matrix to be diagonalizable
		boolean eigenSpaceEqual;
		boolean diagonalizable = true;
		double factor;
		for (int i=0; i<matDim; i++) {
			for (int j=i+1; j<matDim; j++) {
				if ((Math.abs (lambda[i].re - lambda[j].re) < EPSILON) && (Math.abs (lambda[i].im - lambda[j].im) < EPSILON)) {
					// the eigenvalues i and j are equal, so check whether the vectors are equal too
					// we might need a scaling factor, but maybe later [cause they are representatives of the eigenspace]
					factor = V.get0 (0,i).re/V.get0 (0,j).re;
					eigenSpaceEqual = true;
					for (int k=0; k<matDim; k++) {
						eigenSpaceEqual &= (Math.abs (V.get0 (k,i).re - factor*V.get0 (k,j).re) < EPSILON) && (Math.abs (V.get0 (k,i).im - factor*V.get0 (k,j).im) < EPSILON);
					}
					// if the eigenspace is equal they are not diagonalizable
					diagonalizable = diagonalizable && !(eigenSpaceEqual);
				}
			}
		}

		// whatcha say?
		if (!diagonalizable) {
			throw new RuntimeException ("Non-diagonalizable migration matrix encountered.");
		}
		
		//the V matrix should be ok now, get the inverse
		// might have an exception
		try {
			Vinv = Inv.o (V);
		} catch (JampackException e) {
			throw new RuntimeException ("The inverse of V could not be computed.", e);
		}
		// now we should have the inverse
		
		// now we need the VW matrices
		VW = new Zmat[matDim];
		// go through the eigenvalues/vectors and create the corresponding vw
		for (int i=0; i<matDim; i++) {
			// now we create the vw with v being i-th column of V and w being i-th row of Vinv
			Zmat vw = new Zmat (matDim, matDim);
			for (int k=0; k<matDim; k ++) {
				for (int l=0; l<matDim; l++) {
					// hope this works
					vw.put0 (k, l, (new Z()).Times (V.get0 (k, i), Vinv.get0 (i, l)));
				}
			}
			// seem to be good, so store it
			VW[i] = vw;
		}

		// all should be done now
		// check at least whether they has the right length
		assert (lambda.length == matDim);
		assert (VW.length == matDim);		
	}
	
	/// just show something
	public void dump() {
		System.out.println("[VALUES]");
		for (Z lambda : this.lambda) {
			System.out.println(lambda.re + " + i* " + lambda.im );
		}
		System.out.println("[VECTORS - V matrix]");
		for (int i=0; i<this.V.nr; i++) {
			for (int j=0; j<this.V.nc; j++) {
				System.out.print (this.V.get0(i, j).re  + " + i* " +  this.V.get0(i, j).im + "\t");
			}
			System.out.println();
		}
		System.out.println("[inverse]");
		for (int i=0; i<this.Vinv.nr; i++) {
			for (int j=0; j<this.Vinv.nc; j++) {
				System.out.print (this.Vinv.get0(i, j).re  + " + i* " +  this.Vinv.get0(i, j).im + "\t");
			}
			System.out.println();
		}
		System.out.println("[VWs]");
		for (Zmat thisMat : this.VW) {
			System.out.println("{next one}");
			for (int i=0; i<thisMat.nr; i++) {
				for (int j=0; j<thisMat.nc; j++) {
					System.out.print (thisMat.get0(i, j).re  + " + i* " +  thisMat.get0(i, j).im + "\t");
				}
				System.out.println();
			}
		}
	}

	public double[][] getOriginalMatrix () {
		return origMatrix;
	}

	// returns e^(t* M)
	public double[][] getRateMatrixExponential(double t) {
		
		assert (t >= 0d);
		
		// just to make sure
		assert (RateMatrixTools.isProperRateMatrix(this.origMatrix, EPSILON));
		
		// for very short times we get the identity
		if (Math.abs(t) < EPSILON) {
			return Matrix.identity(this.lambda.length, this.lambda.length).getArray();
		}
		
		// for bigger times do the normal exponential
		Zmat matrixSum = new Zmat(this.lambda.length, this.lambda.length);
		// go through eigenstuff
		for (int k=0; k<this.lambda.length; k++) {
			// exponential
			Z exp;
			if (t != Double.POSITIVE_INFINITY) {
				exp = ComplexMath.exponential ((new Z()).Times(t, this.lambda[k]));
			}
			else {
				if (Z.abs (this.lambda[k]) < EPSILON) {
					// eigenvalue zero
					exp = new Z(1d, 0d);
				}
				else {
					// eigenvalue non-zero
					exp = new Z();
				}
			}
			
			// now get the VW on top of matrixSum
			try {
				// hopefully this doesn't mess up
				matrixSum = Plus.o (matrixSum, Times.o (exp, this.VW[k]));
			} catch (JampackException e1) {
				e1.printStackTrace();
				// should not happen
				assert (false);
			}
		}
		// now the matrix should be real
		assert (RateMatrixTools.allEntriesZero (matrixSum.getIm(), EPSILON));
		
		// return the real part
		return matrixSum.getRe();
	}


	@Override
	public int hashCode() {
		// equal orig matrix gives equal hash
		final int prime = 31;
		int result = 1;
		for (double[] row : this.origMatrix) {
			for (double v : row) {
				result = (int)(prime * result + v);
			}
		}
		return result;
	}

	@Override
	public boolean equals(Object obj) {
		// if orig matrices equal, then object is equal
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		MyEigenDecomposition otherEigen = (MyEigenDecomposition) obj;
		if (this.EPSILON != otherEigen.EPSILON) return false;
		if (this.origMatrix.length != otherEigen.origMatrix.length) return false;
		for (int i=0; i<this.origMatrix.length; i++) {
			if (this.origMatrix[i].length != otherEigen.origMatrix[i].length) return false;
			for (int j=0; j<this.origMatrix[i].length; j++) {
				if (this.origMatrix[i][j] != otherEigen.origMatrix[i][j]) return false;
			}
		}
		return true;
	}

	@Override
	public String toString() {
		return "[not implemented yet]";
	}

}
