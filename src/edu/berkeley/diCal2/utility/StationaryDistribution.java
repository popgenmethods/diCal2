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

import java.util.EnumMap;
import java.util.Map.Entry;

import Jama.EigenvalueDecomposition;
import Jama.Matrix;

public class StationaryDistribution {
	public static double[] getStationaryDistribution(double[][] stochasticMatrix) {
		Matrix sMatrix = new Matrix(stochasticMatrix);
		EigenvalueDecomposition eigenDecomp = sMatrix.transpose().eig();
		
		int largestEigenvalueIndex = -1;
		for (int ev = 0; ev < eigenDecomp.getRealEigenvalues().length; ev++) {
			if (largestEigenvalueIndex == -1 || eigenDecomp.getRealEigenvalues()[ev] > eigenDecomp.getRealEigenvalues()[largestEigenvalueIndex]) largestEigenvalueIndex = ev;
		}
		
		Matrix uStationaryDist = eigenDecomp.getV().getMatrix(0, eigenDecomp.getV().getRowDimension() - 1, largestEigenvalueIndex, largestEigenvalueIndex);
		
		double colSum = 0;
		for (int i = 0; i < uStationaryDist.getRowDimension(); i++) colSum += uStationaryDist.get(i, 0);
		
		return uStationaryDist.timesEquals(1 / colSum).getColumnPackedCopy();
	}
	
	
	public static <T extends Enum<T>> EnumMap<T, Double> getStationaryDistribution(EnumMatrix<T, Double> stochasticEnumMatrix) {
		int numAlleles = stochasticEnumMatrix.rowEntrySet().size();
		double[][] sMatrix = new double[numAlleles][numAlleles];
		for (Entry<T, EnumMap<T, Double>> srcEntry : stochasticEnumMatrix.rowEntrySet()) {
			for (Entry<T, Double> dstEntry : srcEntry.getValue().entrySet()) {
				sMatrix[srcEntry.getKey().ordinal()][dstEntry.getKey().ordinal()] = dstEntry.getValue();
			}
		}
		
		double[] sDist = getStationaryDistribution(sMatrix);
		
		EnumMap<T, Double> enumStationaryDist = new EnumMap<T, Double>(stochasticEnumMatrix.getKeyType());
		for (T enumVal : stochasticEnumMatrix.getKeyType().getEnumConstants()) {
			enumStationaryDist.put(enumVal, sDist[enumVal.ordinal()]);
		}
		
		return enumStationaryDist;
	}
}
