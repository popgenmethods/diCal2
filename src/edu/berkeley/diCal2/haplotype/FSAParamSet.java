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

import edu.berkeley.diCal2.csd.EigenParamSet;
import edu.berkeley.diCal2.haplotype.FSAParamSet.FSAUniformMutationParamSet.MutationParameters;
import edu.berkeley.diCal2.utility.CollectionFormat;
import edu.berkeley.diCal2.utility.CopyArray;
import edu.berkeley.diCal2.utility.StationaryDistribution;
import Jama.Matrix;

public class FSAParamSet {
	public interface PSetFactory {
		public FSAUniformMutationParamSet getMutPSet(double theta);
		public FSAUniformRecombinationParamSet getRecPSet(double rho);
	}
	
	public interface PSetFactoryRho {
		public FSAParamSet getPSet(double rho);
	}
	
	/// This has the same mutation rate and mutation matrix at every locus
	public static class FSAUniformMutationParamSet {
		public final static Matrix PIM_MATRIX_2 = getPIMMatrix(2);
		
		private final int numLoci;
		private final int numAlleles;

		public static class MutationParameters {
			private final double mutRate;
			private final Matrix mutMatrix;
			
			// stationary probabilities (derived)
			private final double[] stationaryProbability;
			
			// PIM description (derived)
			private final Matrix npimMutMatrix;
			private final double[] pimMutArray;
			private final double pimProbability;
			
			public MutationParameters(double mutRate, Matrix mutMatrix) {
				this.mutRate = mutRate;
				this.mutMatrix = mutMatrix;
				
				int numAlleles = mutMatrix.getRowDimension();
				
				this.pimMutArray = new double[numAlleles];
				this.npimMutMatrix = new Matrix(numAlleles, numAlleles);
				double pimProbability = 0;
				
				// decompose the mutation matrix
				for (int dstAllele = 0; dstAllele < numAlleles; dstAllele++) {
					double minSrcProb = 1;
					for (int srcAllele = 0; srcAllele < numAlleles; srcAllele++) {
						minSrcProb = Math.min(minSrcProb, mutMatrix.get(srcAllele, dstAllele));
					}

					this.pimMutArray[dstAllele] = minSrcProb;
					pimProbability += minSrcProb;
				}

				// re-normalize
				for (int dstAllele = 0; dstAllele < numAlleles; dstAllele++) {
					for (int srcAllele = 0; srcAllele < numAlleles; srcAllele++) {
						this.npimMutMatrix.set(srcAllele, dstAllele, (mutMatrix.get(srcAllele, dstAllele) - pimMutArray[dstAllele]) / (1 - pimProbability));
					}

					this.pimMutArray[dstAllele] /= pimProbability;
				}
				
				this.pimProbability = pimProbability;
				
				stationaryProbability = StationaryDistribution.getStationaryDistribution(mutMatrix.getArray());
			}
			
			public double getMutationRate() { return mutRate; }
			public Matrix getMutationMatrix() { return mutMatrix; }

			// PIM description
			public double getPIMProbability() { return pimProbability; }
			public double[] getPIMMutationArray() { return pimMutArray; }
			public Matrix getNonPIMMutationMatrix() { return npimMutMatrix; }

			// stationary probability
			public double getStationaryProbability(int allele) { return stationaryProbability[allele]; }
			
			public String toString() {
				return Double.toString(mutRate);
			}
		}
		
//		private final MutationParameters[] mutationParameters;
		private final MutationParameters uniformMutationParameters;
		
//		FSAUniformMutationParamSet(int numLoci, int numAlleles, double mutRate, Matrix mutMatrix) {
//			this.numLoci = numLoci;
//			this.numAlleles = numAlleles;
//			this.uniformMutationParameters = new MutationParameters[numLoci];
//			
//			for (int locus = 0; locus < numLoci; locus++) {
//				this.uniformMutationParameters[locus] = new MutationParameters(mutRate[locus], mutMatrix[locus]);
//			}
//		}
		
		public FSAUniformMutationParamSet(int numLoci, int numAlleles, double mutRate, Matrix mutMatrix) {
			this.numLoci = numLoci;
			this.numAlleles = numAlleles;
			
//			MutationParameters sMutationParameters = new MutationParameters(mutRate, mutMatrix);

			this.uniformMutationParameters = new MutationParameters(mutRate, mutMatrix);
		}
		
		public static Matrix getPIMMatrix(int numAlleles) {
			double[][] pimRawMatrix = new double[numAlleles][numAlleles];
			for (int as = 0; as < numAlleles; as++) {
				for (int ad = 0; ad < numAlleles; ad++) {
					pimRawMatrix[as][ad] = 1.0 / numAlleles;
				}
			}
			
			return new Matrix(pimRawMatrix);
		}

		public int numLoci() { return numLoci; }
		public int numAlleles() { return numAlleles; }
		
		// all same at every locus
		public MutationParameters getMutationParameters(int locus) { return this.uniformMutationParameters; }
		
		// convenience methods
		// all same at every locus
		public double getMutationRate(int locus) { return uniformMutationParameters.getMutationRate(); }
		public Matrix getMutationMatrix(int locus) { return uniformMutationParameters.getMutationMatrix(); }
		
		// PIM description
		// all same at every locus
		public double getPIMProbability(int locus) { return uniformMutationParameters.getPIMProbability(); }
		public double[] getPIMMutationArray(int locus) { return uniformMutationParameters.getPIMMutationArray(); }
		public Matrix getNonPIMMutationMatrix(int locus) { return uniformMutationParameters.getNonPIMMutationMatrix(); }
		
		// stationary probability
		// all same at every locus
		public double getStationaryProbability(int locus, int allele) { return uniformMutationParameters.getStationaryProbability(allele); }
		
		public String toString() {
//			return CollectionFormat.formatArray(this.uniformMutationParameters, ",", "[", "]");
			return "[" + this.uniformMutationParameters.toString() + "]";
		}
	}
	
	public static class FSAUniformRecombinationParamSet {
		private final int numLoci;
//		private final double[] recRate;
		private final double recRate;
		
//		public FSAUniformRecombinationParamSet(int numLoci, double[] recRate) {
//			this.numLoci = numLoci;
//			this.recRate = recRate;
//		}
		
		public FSAUniformRecombinationParamSet(int numLoci, double recRate) {
			this.numLoci = numLoci;
			this.recRate = recRate;
		}
		
		// unfortunately, we need to suspend this for now
//		public static FSAUniformRecombinationParamSet getReducedParamSet(LocusSet sites, double[] originalRecRate) {
//			double[] adjustedRecRate = new double[sites.size() - 1];
//			
//			int prevOldLocus = -1;
//			int newIdx = 0;
//
//			for (int locus : sites) {
//				if (prevOldLocus != -1) {
//					int totalRate = 0;
//					for (int l = prevOldLocus; l < locus; l++) {
//						totalRate += originalRecRate[l];
//					}
//				
//					adjustedRecRate[newIdx++] = totalRate;
//				}
//				
//				prevOldLocus = locus;
//			}
//			
//			return new FSAUniformRecombinationParamSet(sites.size(), adjustedRecRate);
//		}
		
		// unfortunately, we need to suspend this for now
//		public static FSAUniformRecombinationParamSet getReducedParamSet(LocusSet sites, double originalRecRate) {
//			return getReducedParamSet(sites, CopyArray.nCopies(sites.size() - 1, originalRecRate));
//		}
		
		public int numLoci() { return numLoci; }
		
		public double getRecombinationRate(int startLocus, int endLocus) {
			assert (startLocus <= endLocus);
			// just linear
			return (endLocus - startLocus) * this.recRate;
			
//			if (endLocus - startLocus == 1) {
//				return recRate;
//			} else {
//				double totalRate = 0;
//				for (int l = startLocus; l < endLocus; l++) {
//					totalRate += recRate[l];
//				}
//				
//				return totalRate;
//			}
		}
		
		public String toString() {
//			return CollectionFormat.formatArray(recRate, ",", "[", "]");
			return "[" + this.recRate + "]";
		}
	}
	
	private final FSAUniformMutationParamSet mutPSet;
	private final FSAUniformRecombinationParamSet recPSet;
	
	public FSAParamSet(FSAUniformMutationParamSet mutPSet, FSAUniformRecombinationParamSet recPSet) {
		this.mutPSet = mutPSet;
		this.recPSet = recPSet;
	}
	
//	public FSAParamSet(int numLoci, int numAlleles, double[] mutRate, Matrix[] mutMatrix, double[] recRate) {
//		this(new FSAUniformMutationParamSet(numLoci, numAlleles, mutRate, mutMatrix), new FSARecombinationParamSet(numLoci, recRate));
//	}
	
	public FSAParamSet(int numLoci, int numAlleles, double mutRate, Matrix mutMatrix, double recRate) {
		this(new FSAUniformMutationParamSet(numLoci, numAlleles, mutRate, mutMatrix), new FSAUniformRecombinationParamSet(numLoci, recRate));
	}

	public FSAUniformMutationParamSet getMutationParamSet() { return mutPSet; }
	public FSAUniformRecombinationParamSet getRecombinationParamSet() { return recPSet; }

	public int numLoci() { return mutPSet.numLoci(); }
	public int numAlleles() { return mutPSet.numAlleles(); }

	public MutationParameters getMutationParameters(int locus) { return mutPSet.getMutationParameters(locus); }

	public double getMutationRate(int locus) { return mutPSet.getMutationRate(locus); }
	public Matrix getMutationMatrix(int locus) { return mutPSet.getMutationMatrix(locus); }

	// PIM description
	public double getPIMProbability(int locus) { return mutPSet.getPIMProbability(locus); }
	public double[] getPIMMutationArray(int locus) { return mutPSet.getPIMMutationArray(locus); }
	public Matrix getNonPIMMutationMatrix(int locus) { return mutPSet.getNonPIMMutationMatrix(locus); }

	// stationary probability
	public double getStationaryProbability(int locus, int allele) { return mutPSet.getStationaryProbability(locus, allele); }
	
	// the total recombination rate between two loci
	public double getRecombinationRate(int startLocus, int endLocus) { return recPSet.getRecombinationRate(startLocus, endLocus); }
	
	public String toString() {
		return "mutation: " + mutPSet + ", recombination: " + recPSet;
	}
}
