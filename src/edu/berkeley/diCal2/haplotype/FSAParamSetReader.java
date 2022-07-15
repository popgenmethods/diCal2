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

package edu.berkeley.diCal2.haplotype;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.Reader;
import java.util.ArrayList;

import edu.berkeley.diCal2.utility.RateMatrixTools;
import Jama.Matrix;

// This class reads in an FSA Param set from a file
public class FSAParamSetReader {

	public static class OneLocusParams {
		public final double theta;
		public final double rho;
		public final Matrix mutMatrix;
		public final int numAlleles;
		
		public OneLocusParams(Reader r, double EPSILON, boolean renormalizeMutationMatrix) throws IOException {
			BufferedReader bReader = new BufferedReader(r);
			
			
			ArrayList<String> lines = new ArrayList<String>();
			String line;
			while ((line = bReader.readLine()) != null) {
				if (line.trim().startsWith("#") || line.trim().equals("")) continue;
				lines.add(line);
			}
			bReader.close();
			
			double thetaVal = Double.parseDouble(lines.get(0));
			lines.remove(0);
			double rhoVal = Double.parseDouble(lines.get(0));
			lines.remove(0);
			
			int numAlleles = lines.size();
			if (numAlleles <= 1) throw new IOException("Improperly formatted mutation matrix.");
			
			
			double[][] mutationMatrix = new double[numAlleles][numAlleles];
			for (int i = 0; i < numAlleles; i++) {
				String[] currLine = lines.get(i).split("(\t| )+");
				if (currLine.length != numAlleles) throw new IOException("Improperly formatted mutation matrix.");
				
				double rowSum = 0d;
				for (int j = 0; j < numAlleles; j++) {
					rowSum += mutationMatrix[i][j] = Double.parseDouble(currLine[j]);
				}
				if (renormalizeMutationMatrix) {
					for (int j = 0; j < numAlleles; j++) {
						mutationMatrix[i][j] /= rowSum;
					}
				}
			}
			
			if (!RateMatrixTools.isProperStochaticMatrix(mutationMatrix, EPSILON)) throw new IOException("Mutation matrix is not a stochastic matrix.");

			this.numAlleles = numAlleles;
			this.theta = thetaVal;
			this.rho = rhoVal;
			this.mutMatrix = new Matrix(mutationMatrix);
		}
	}
	
	public static FSAParamSet readParams(Reader r, int numLoci, double EPSILON) throws IOException {
		return readParams(r, numLoci, EPSILON, false);
	}
	
	public static FSAParamSet readParams(Reader r, int numLoci, double EPSILON, boolean renormalizeMutationMatrix) throws IOException {
		OneLocusParams par = new OneLocusParams(r, EPSILON, renormalizeMutationMatrix);
		return new FSAParamSet(numLoci, par.numAlleles, par.theta, par.mutMatrix, par.rho);
	}

}
