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
import java.io.IOException;
import java.io.Reader;

import edu.berkeley.diCal2.csd.auxiliary.MyEigenDecomposition;
import edu.berkeley.diCal2.haplotype.FSAParamSet;
import edu.berkeley.diCal2.haplotype.FSAParamSetReader.OneLocusParams;
import Jama.Matrix;


public class EigenParamSet extends FSAParamSet {

	private final MyEigenDecomposition mutMatrixEigenDecomp;
	
	public EigenParamSet(int numLoci, int numAlleles, double mutRate,
			Matrix mutMatrix, double recRate) {
		super(numLoci, numAlleles, mutRate, mutMatrix, recRate);
		
		this.mutMatrixEigenDecomp = new MyEigenDecomposition(mutMatrix.getArray(), UberDemographyCore.ONELOCUS_EPSILON);		
	}

	public static EigenParamSet ReadEigenParamSet(Reader r, int numLoci, double EPSILON, boolean renormalizeMutationMatrix) throws IOException {
		OneLocusParams par = new OneLocusParams(r, EPSILON, renormalizeMutationMatrix);
		return new EigenParamSet(numLoci, par.numAlleles, par.theta, par.mutMatrix, par.rho);
	}

	public static EigenParamSet ReadEigenParamSet(Reader r, int numLoci, double EPSILON) throws IOException {
		return ReadEigenParamSet(r, numLoci, EPSILON, false);
	}

	public MyEigenDecomposition getMutMatrixEigenDecomp() {
		return mutMatrixEigenDecomp;
	}
}
