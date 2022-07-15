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
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;
import java.util.concurrent.ConcurrentHashMap;

import edu.berkeley.diCal2.csd.DemoState.DemoStateCollection;
import edu.berkeley.diCal2.csd.DemoState.DemoStateFactory;
import edu.berkeley.diCal2.csd.TrunkProcess.TrunkProcessFactory;
import edu.berkeley.diCal2.csd.auxiliary.MyEigenDecomposition;
import edu.berkeley.diCal2.demography.Demography;
import edu.berkeley.diCal2.haplotype.FSAHaplotype;
import edu.berkeley.diCal2.haplotype.FSAParamSet;
import edu.berkeley.diCal2.haplotype.FSAXHaplotype;
import edu.berkeley.diCal2.haplotype.FSAXTypeConfig.FSAReference;
import edu.berkeley.diCal2.haplotype.GeneticType.GeneticTypeMultiplicity;
import edu.berkeley.diCal2.haplotype.HFSAXFullHaplotype.HFSAXFullHaplotypeShell;
import edu.berkeley.diCal2.haplotype.ReadSequences;
import edu.berkeley.diCal2.maximum_likelihood.MetaOptimization;
import edu.berkeley.diCal2.utility.LogSum;
import edu.berkeley.diCal2.utility.MemSaveArrays.FiveDimShortArray;
import edu.berkeley.diCal2.utility.SimplePair;
import edu.berkeley.diCal2.utility.StationaryDistribution;
import gnu.trove.set.TIntSet;
import gnu.trove.set.hash.TIntHashSet;

public abstract class HmmStepHandler<H extends FSAHaplotype> {

	protected final HashMap<H, Double> stationaryMap = new HashMap<H, Double>();
	protected final EigenParamSet pSet;

	
	public HmmStepHandler (EigenParamSet pSet) {
		this.pSet = pSet;
	}
	
	// decoding
	public abstract int numHmmSteps ();
	
	// used to check if we have the emission table right
	public abstract int numNonEmissionLoci();
	
	public abstract int getTransitionIndex (int decodeLocusIdx, boolean reverseDirection);
	
	public abstract double[] getMutationRate (int decodeLocusIdx); // assumes all loci within an hmmStep have the same mutation rate
	
	public abstract SMCSDTransitionList computeTransitions (UberDemographyCore core, double recoScaling);
	
	public abstract double getLogEmissionProb (double[][] singleBaseLogEmission, int hmmStep, H trunkHap, int trunkPresentDeme, H addHap);
	public abstract void updateEmissionTable (double condProb, double[][] emissionTable, int hmmStep, H trunkHap, H addHap, double[][] logEmissionArray);
	
	public void updateTheta(double[] theta){
		//Does nothing, but some classes override this that need it
	}
	
	public double computeLogStationaryProbability (H additionalHaplotype) {
		// did we already compute it?
		Double value = stationaryMap.get(additionalHaplotype);
		if (value == null) {
			// apparently not, so compute it

			// but everything in log scale
			value = 0d;
			// go though loci
			for (int l=0; l<this.pSet.numLoci(); l++) {
				int allele = additionalHaplotype.getAllele(l);
				// for missing allele we would add 0 [log(1)]
				if (allele != ReadSequences.MISSING_ALLELE) {
					value += Math.log (this.pSet.getStationaryProbability(l, allele));
				}
			}
			// should be it
			
			// put it in the map
			stationaryMap.put(additionalHaplotype, value);
		}
		
		// return it
		return value;
	}
	
	abstract public boolean isValidHaplotype (H additionalHaplotype);
	abstract public boolean isValidConfig (HaplotypeConfiguration<H> config);
	
	// just show some things
	abstract public void printInfo (PrintStream outStream);	
	
	protected static <H extends FSAHaplotype> double easyGetEmissionProb (double[][] singleBaseEmission, int locus, H trunkHap, H addHap) {
		return singleBaseEmission[trunkHap.getAllele(locus)][addHap.getAllele(locus)];
	}
	
	protected static <H extends FSAHaplotype> void easyUpdateEmissionTable (double condProb, double[][] singleBaseEmission, int locus, H trunkHap, H addHap) {
		singleBaseEmission[trunkHap.getAllele(locus)][addHap.getAllele(locus)] += condProb;
	}
	
	protected static <H extends FSAHaplotype> SMCSDTransitionList easyComputeTransitions(UberDemographyCore core, List<Double> recomRates, double recoScaling) {
		
		// should never happen
		assert (!core.isConditionalConfigEmpty());
		
		// storage or the transitions
		SMCSDTransitionList transitionList = new SMCSDTransitionList ();
		
		// go through all rates present and compute
		for (int i = 0; i < recomRates.size(); i++){
			
			UberDemographyCore.RecoLogProbs recomLogProbs = core.getRecoLogProbs(recomRates.get(i) * recoScaling);
			
			transitionList.setLogNoReco (i, recomLogProbs.logNoReco);
			transitionList.setLogReco (i, recomLogProbs.logReco);
		}
		
		// return list object
		return transitionList;
	}

	
	public interface HasPartiallyMissing {
		public int numPartiallyMissing();
	}
	
	public static abstract class HmmLocusStepHandler<H extends FSAHaplotype> extends HmmStepHandler<H> {

		public HmmLocusStepHandler(EigenParamSet pSet) {
			super(pSet);
		}

		public abstract int getLocus(int decodeLocusIdx);
	}
	
	public static class SingleStepStupidMissingAlleles<H extends FSAHaplotype> extends SingleStepLocusTransitionMap<H> implements HasPartiallyMissing{

		public SingleStepStupidMissingAlleles(EigenParamSet pSet) {
			super(pSet);
		}

		@Override
		public double getLogEmissionProb(double[][] singleBaseEmission, int hmmStep, H trunkHap, int trunkPresentDeme, H addHap) {
			
			//addHap allele is missing, so do nothing
			if(addHap.getAllele(hmmStep) == ReadSequences.MISSING_ALLELE){
				return 0;
			}
			
			
			//trunkHap allele is missing, use stationary distribution
			if(trunkHap.getAllele(hmmStep) == ReadSequences.MISSING_ALLELE){
				return Math.log(this.pSet.getStationaryProbability(hmmStep, addHap.getAllele(hmmStep)));
			}
			
			return super.getLogEmissionProb(singleBaseEmission, hmmStep, trunkHap, trunkPresentDeme, addHap);
		}


		@Override
		public void updateEmissionTable(double condProb,
				double[][] emissionTable, int hmmStep, H trunkHap, H addHap, double[][] logEmissionProbs) {
			
			// addHap allele is missing, so do nothing
			if(addHap.getAllele(hmmStep) == ReadSequences.MISSING_ALLELE){
				return;
			}

			if (trunkHap.getAllele(hmmStep) == ReadSequences.MISSING_ALLELE) {
				
				int locus = hmmStep;
				int numAlleles = emissionTable.length;
				
				double[] logStationary = new double[numAlleles];
				for (int i = 0; i < numAlleles; i++) logStationary[i] = Math.log(this.pSet.getStationaryProbability(locus, i));
				
				// conditional probability of hidden state and trunk allele, given the data
				double[] alleleCondProbs = getStupidStationaryTrunkAllelePosterior(addHap.getAllele(locus), condProb, logStationary, logEmissionProbs);
				
				for (int trunkAllele = 0; trunkAllele < numAlleles; trunkAllele++) {
					emissionTable[trunkAllele][addHap.getAllele(locus)] += alleleCondProbs[trunkAllele];
				}
			} else {
				super.updateEmissionTable(condProb, emissionTable, hmmStep, trunkHap, addHap, null);
			}
		}
		
		@Override
		public int numNonEmissionLoci() {
			// what to do here?
			assert false;
			return 0;
		}

		@Override
		public int numPartiallyMissing() {
			// what to do here?
			assert false;
			return 0;
		}

	}
	
	public static class SingleStepClevererImpute<H extends FSAHaplotype> extends SingleStepMissingAllelesTransitionMap<H> {

		final DemoConfiguration<H> nonmissingHaps;
		
		public SingleStepClevererImpute(
				EigenParamSet pSet,
				DemoConfiguration<H> configs,
				CoreCache<H> coreCache) {
			super(pSet, configs, coreCache);
			
			this.nonmissingHaps = new DemoConfiguration<H>(configs.getNumberDemes(), configs.numLoci(), configs.numAlleles());
			for (int deme = 0; deme < configs.getNumberDemes(); deme++) {
				for (GeneticTypeMultiplicity<H> hapCount : configs.getPopulation(deme)) {
					H hap = hapCount.geneticType;
					boolean hasMissing = false;
					for (int locus = 0; locus < hap.getNumLoci(); locus++) {
						hasMissing |= (hap.getAllele(locus) == ReadSequences.MISSING_ALLELE);
					}
					if (!hasMissing) {
						this.nonmissingHaps.getPopulation(deme).adjustType(hap, hapCount.multiplicity);
					}
				}
			}
		}

		@Override
		protected double[] getImputedAlleles(double[][] singleBaseEmission,
				int hmmStep, H trunkGuy, int trunkPresentDeme, H addHap) {
			if (trunkGuy.getAllele(hmmStep) != ReadSequences.MISSING_ALLELE) {
				// returns delta function in this case
				return super.getImputedAlleles(singleBaseEmission, hmmStep, trunkGuy,
				trunkPresentDeme, addHap);
			}
			
			if (this.nonmissingHaps.getTotalGeneticTypes() == 0) {
				return super.getImputedAlleles(singleBaseEmission, hmmStep, trunkGuy, trunkPresentDeme, addHap);
			}
			
			UberDemographyCore nonmissingCore = coreCache.getCore(this.nonmissingHaps, trunkPresentDeme);
			SMCSDemo<H> smcsd = new SMCSDemo<H>(this.nonmissingHaps, nonmissingCore, pSet, 1d, new SingleStepStupidMissingAlleles<H>(pSet));
			smcsd.fillForwardBackward(trunkGuy);
			
			double[] toReturn = new double[this.pSet.numAlleles()];
			
			// [demoState][fromAllele][toAllele]
			double[][][] demoStateLogEmissions = smcsd.core.getLogEmission(new double[] {pSet.getMutationRate(hmmStep)}, pSet.getMutMatrixEigenDecomp());
			
			for (int demoState=0; demoState<smcsd.core.numDemoStates(); demoState++) {
				int currDemeIdx = smcsd.core.getDemoStateCollection().getPresentDeme(demoState);
				HaplotypeConfiguration<H> currPop = smcsd.getPresentConfig().getPopulation(currDemeIdx);
				
				double[] currDemoStateAllelePosteriors = new double[this.pSet.numAlleles()];
				
				int hapIndex = 0;
				for (GeneticTypeMultiplicity<H> nonmissingHap : currPop) {
					
					double nonmissingHapStateProb = smcsd.getPosteriorProb(hmmStep, demoState, hapIndex);
					int nonmissingAllele = nonmissingHap.geneticType.getAllele(hmmStep);
					
					currDemoStateAllelePosteriors[nonmissingAllele] += nonmissingHapStateProb;
					
					// increase index
					hapIndex++;
				}
				
				for (int nonmissingAllele = 0; nonmissingAllele < this.pSet.numAlleles(); nonmissingAllele++) {
					for (int missingTrunkAllele = 0; missingTrunkAllele < this.pSet.numAlleles(); missingTrunkAllele++) {						
						toReturn[missingTrunkAllele] += currDemoStateAllelePosteriors[nonmissingAllele] * Math.exp(demoStateLogEmissions[demoState][nonmissingAllele][missingTrunkAllele]);
					}
				}
			}
			
			return toReturn;
		}
	}
	
	// Takes care of missing alleles
	public static class SingleStepMissingAllelesTransitionMap<H extends FSAHaplotype> extends SingleStepLocusTransitionMap<H> implements HasPartiallyMissing{
		
		private final DemoConfiguration<H> configs;
		protected final CoreCache<H> coreCache;
		
		public static class CoreCache <H extends FSAHaplotype> {
			private final ConcurrentHashMap<SimplePair<List<Integer>,Integer>, UberDemographyCore> cache;
			private final Demography demo;
			private final TrunkProcessFactory trunkFactory;
			private final DemoStateFactory demoStateFactory;
			private final double renormalizeEpsilon;
			private final boolean useEigenCore;
			
			public CoreCache(Demography demo, TrunkProcessFactory trunkFactory, DemoStateFactory demoStateFactory, double renormalizeEpsilon, boolean useEigenCore) {
				super();
				this.demo = demo;
				this.cache = new ConcurrentHashMap<SimplePair<List<Integer>,Integer>, UberDemographyCore>();
				this.trunkFactory = trunkFactory;
				this.demoStateFactory = demoStateFactory;
				this.renormalizeEpsilon = renormalizeEpsilon;
				
				this.useEigenCore = useEigenCore;
			}

			public UberDemographyCore getCore (DemoConfiguration<H> presentConfig, int additionalDeme) {
				SimplePair<List<Integer>, Integer> hashKey= getCoreKey(presentConfig.getSampleSizes(), additionalDeme);
				
				if (!cache.containsKey(hashKey)) {
					DemoStateCollection demoStates = demoStateFactory.getDemoStates(demo, presentConfig, additionalDeme);
					TrunkProcess trunk = trunkFactory.getTrunk(presentConfig.getSampleSizes(), demoStates);
					cache.putIfAbsent(hashKey, UberDemographyCore.getDemographyCore(trunk, demoStates, additionalDeme, true, renormalizeEpsilon, this.useEigenCore));
				}
				
				return cache.get(hashKey);
			}
			
			private static SimplePair<List<Integer>, Integer> getCoreKey(int[] sampleSizes, int addDeme){
				List<Integer> sampleSizeList = new ArrayList<Integer>();
				for (int size : sampleSizes) sampleSizeList.add(size);
				return new SimplePair<List<Integer>, Integer>(sampleSizeList, addDeme);
			}

			public Demography getDemo() {
				return this.demo;
			}

			public boolean isEigenCore() {
				return this.useEigenCore;
			}
		}

		
		public SingleStepMissingAllelesTransitionMap(EigenParamSet pSet, DemoConfiguration<H> configs, CoreCache<H> coreCache) {
			super(pSet);
			this.configs= configs;
			this.coreCache= coreCache;
		}
		

		public static <H extends FSAHaplotype> double[] getTrunkHapAlleleProbs(int locus, H trunkHap, int trunkPresentDeme,
				DemoConfiguration<H> configs, CoreCache<H> coreCache, EigenParamSet pSet) {
			if (trunkHap.getAllele(locus) != ReadSequences.MISSING_ALLELE){
				double[] toReturn = new double[pSet.numAlleles()];
				toReturn[trunkHap.getAllele(locus)] = 1d;
				return toReturn;
			}
			
			int numDemes= configs.getNumberDemes();
			int[][] alleleCounts= new int[numDemes][pSet.numAlleles()];
			int[] totalCount= new int[numDemes];
			DemoConfiguration<H> nonMissingAtLocus = new DemoConfiguration<H>(numDemes, 2, pSet.numAlleles());
			
			for(int deme= 0; deme < numDemes; deme++){
				for(GeneticTypeMultiplicity<H> type : configs.getPopulation(deme)){
					if(type.geneticType.getAllele(locus) != ReadSequences.MISSING_ALLELE){
						alleleCounts[deme][type.geneticType.getAllele(locus)] += type.multiplicity;
						totalCount[deme] += type.multiplicity;
						nonMissingAtLocus.adjustHaplotypeMultiplicity(deme, type.geneticType, type.multiplicity);
					}
				}
			}
			if(nonMissingAtLocus.getTotalGeneticTypes() == 0){
				double[] toReturn = new double[pSet.numAlleles()];
				for (int allele = 0; allele < pSet.numAlleles(); allele++) {
					toReturn[allele] = pSet.getStationaryProbability(locus, allele);
				}
				return toReturn;
			}
			UberDemographyCore core= coreCache.getCore(nonMissingAtLocus, trunkPresentDeme);
			
			return core.imputeLocus(alleleCounts, core.getLogEmission(new double[] {pSet.getMutationRate(locus)}, pSet.getMutMatrixEigenDecomp()));
		}
		
		public static <H extends FSAHaplotype> double oneLocusEmissionProb(double[][] singleBaseEmission,
				int locus, double[] imputedAllele, int trunkPresentDeme, H addHap) {
			// addHap allele is missing, so do nothing
			if(addHap.getAllele(locus) == ReadSequences.MISSING_ALLELE){
				return 0;
			}
			
			//trunkHap allele is missing, sample uniformly from other haps
			double totalImputedProb = 0;
			double prob= 0;
			for(int fromAllele= 0; fromAllele < singleBaseEmission.length; fromAllele++){
				prob += Math.exp(singleBaseEmission[fromAllele][addHap.getAllele(locus)]) * imputedAllele[fromAllele];
				totalImputedProb += imputedAllele[fromAllele];
			}
			assert Math.abs(totalImputedProb - 1d) < UberDemographyCore.ONELOCUS_EPSILON;
			return Math.log(prob);
		}
		
		@Override
		public final double getLogEmissionProb(double[][] singleBaseEmission,
				int hmmStep, H trunkHap, int trunkPresentDeme, H addHap) {
			double[] imputedAlleles = getImputedAlleles(singleBaseEmission, hmmStep, trunkHap, trunkPresentDeme, addHap);
			return oneLocusEmissionProb(singleBaseEmission, hmmStep, imputedAlleles, trunkPresentDeme, addHap);
		}

		protected double[] getImputedAlleles(double[][] singleBaseEmission,
				int hmmStep, H trunkHap, int trunkPresentDeme, H addHap) {
			return getTrunkHapAlleleProbs(hmmStep, trunkHap, trunkPresentDeme, configs, coreCache, pSet);
		}

		@Override
		public void updateEmissionTable(double condProb,
				double[][] emissionTable, int hmmStep, H trunkHap, H addHap,
				double[][] logEmissionArray) {
			throw new RuntimeException("EM not implemented for missing alleles step handler. (use the stupid missing alleles handler instead)");
		}
		
		@Override
		public int numNonEmissionLoci() {
			// what to do here?
			assert false;
			return 0;
		}


		@Override
		public int numPartiallyMissing() {
			// what to do here?
			assert false;
			return 0;
		}

		
		
	}
	
	public static class SingleStepLocusTransitionMap<H extends FSAHaplotype> extends HmmLocusStepHandler<H> {		

		public SingleStepLocusTransitionMap (EigenParamSet pSet) {
			
			super(pSet);
			
			// create storage for rates nd map
			this.recomRates = new ArrayList<Double>();
			this.locusToRecomRateIdxMap = new int[pSet.numLoci()-1];
			
			// go through loci
			for (int locus = 0; locus < pSet.numLoci()-1; locus++){
				
				// what rate do we have here?
				double currRate = this.pSet.getRecombinationRate(locus, locus+1);
				
				// is the rate already in storage?
				boolean foundRate = false;				
				for (int i = 0; i < this.recomRates.size(); i++){
					double listRate = this.recomRates.get(i);
					if (Math.abs ((currRate - listRate)/listRate) < UberDemographyCore.RELATIVE_ERROR) {
						// yes it is, remember where
						this.locusToRecomRateIdxMap[locus] = i;
						foundRate = true;
						break;
					}
				}
				
				// we went through and its not in there, so put it in there
				if (!foundRate){
					this.locusToRecomRateIdxMap[locus] = this.recomRates.size();
					this.recomRates.add(currRate);
				}
			}
			// this should be done
		}
		
		@Override
		public int getTransitionIndex(int decodeLocusIdx, boolean reverseDirection) {
			assert (!reverseDirection || (decodeLocusIdx > 0));
			assert (reverseDirection || (decodeLocusIdx < this.locusToRecomRateIdxMap.length));
			return this.locusToRecomRateIdxMap[decodeLocusIdx - (reverseDirection ? 1 : 0)];
		}

		@Override
		public int numNonEmissionLoci() {
			// no non-emission here
			return 0;
		}

		@Override
		public SMCSDTransitionList computeTransitions(UberDemographyCore core, double recoScaling) {
			return easyComputeTransitions (core, this.recomRates, recoScaling);
		}

		@Override
		public int numHmmSteps() {
			return this.locusToRecomRateIdxMap.length+1;
		}

		@Override
		public int getLocus(int decodeLocusIdx) {
			return decodeLocusIdx;
		}



		@Override
		public double[] getMutationRate(int decodeLocusIdx) {
			double[] toReturn= new double[1];
			toReturn[0]= pSet.getMutationRate(decodeLocusIdx);
			return toReturn;
		}
		
		@Override
		public double getLogEmissionProb(double[][] singleBaseEmission, int hmmStep,
				H trunkHap, int trunkPresentDeme, H addHap) {
			return easyGetEmissionProb(singleBaseEmission, hmmStep, trunkHap, addHap);
		}
		
		@Override
		public void updateEmissionTable(double condProb, double[][] emissionTable,
				int hmmStep, H trunkHap, H addHap, double[][] logEmissionArray) {
			easyUpdateEmissionTable(condProb, emissionTable, hmmStep, trunkHap, addHap);
		}
		
		@Override
		public boolean isValidHaplotype(H additionalHaplotype) {
			// always valid
			return true;
		}
		
		@Override
		public boolean isValidConfig (HaplotypeConfiguration<H> config) {
			// evrytihng's valid
			return true;
		}

		@Override
		public void printInfo (PrintStream outStream) {
			System.out.println ("# Total number loci: " + this.pSet.numLoci());
		}

		private final int[] locusToRecomRateIdxMap;
		private final List<Double> recomRates;
		
	}
	
	public static class SimpleLocusSkipper<H extends FSAHaplotype> extends HmmLocusStepHandler<H> {
		
		// constructor
		public SimpleLocusSkipper (List<H> validSequenceList, EigenParamSet pSet, TreeSet<Integer> additionalDecodeLoci) {	//pSet contains locus varying recombination
			super (pSet);
			
			// check whether all recombination rates, all mutation rates and all mutation matrices in pSet are equal
			// and whether the mutation matrix is symmetric
			this.numLoci = pSet.numLoci();
			assert (this.numLoci > 1);
			this.numAlleles = pSet.numAlleles();
			
			this.recoRate = pSet.getRecombinationRate(0, 1);
			this.mutRate= new double[1];
			this.mutRate[0] = pSet.getMutationRate(0);
			double[][] mutMatrix = pSet.getMutationMatrix(0).getArray();
			this.firstEigenMutMatrix = new MyEigenDecomposition (mutMatrix, UberDemographyCore.ONELOCUS_EPSILON);
			assert(isUniformMutationMatrix(mutMatrix, UberDemographyCore.ONELOCUS_EPSILON));
			// the stationary distribution for the mutation should be the uniform distribution
			assert (isStationaryUniformDistribution (pSet, 0, UberDemographyCore.ONELOCUS_EPSILON));
				
			// assert some stuff
//			assert (SMCSDemo.matrixEqual (mutMatrix, pSet.getMutationMatrix(0).getArray(), UberDemographyCore.ONELOCUS_EPSILON, UberDemographyCore.RELATIVE_ERROR));
			// with the new type for all parameters same everywhere, this is given
//			assert (SMCSDemo.checkMutationMatricesEqual(pSet, UberDemographyCore.ONELOCUS_EPSILON, UberDemographyCore.RELATIVE_ERROR));
//			assert (checkMutationAndRecombinationRatesAllEqual(pSet, UberDemographyCore.ONELOCUS_EPSILON));
			
			// get the total population
			this.validSequences = new HashSet<H>(validSequenceList);
			
			
			// find the segregating sites (will become the decode loci [more or less])
			TreeSet<Integer> decodeLociSet = getSegregatingSites (this.validSequences);
			this.numSegSites = decodeLociSet.size();
					
			// for convenience we add 0 and the last locus
			decodeLociSet.add(0);
			decodeLociSet.add(pSet.numLoci() - 1);
			
			// add additional decode loci
			if (additionalDecodeLoci != null) decodeLociSet.addAll(additionalDecodeLoci);
			
			// now find the actual decode loci
			this.decodeLoci = new ArrayList<Integer>();

			// iterate over segregating sites
			Iterator<Integer> itSegSites = decodeLociSet.iterator();
			// first one should be zero
			int rightSegSite = itSegSites.next();
			int leftSegSite = rightSegSite;
			assert (rightSegSite == 0);
			
			// since our following method is right inclusive, we ad this guy here
			decodeLoci.add (0);
			
			// loop through it
			while (itSegSites.hasNext()) {
				// update 
				leftSegSite = rightSegSite;
				rightSegSite = itSegSites.next();
				
				assert (leftSegSite < rightSegSite);
				// do stuff
				decodeLoci.addAll (getPowerOfTwoStepsBetween  (leftSegSite, rightSegSite));
			}
			assert(rightSegSite == pSet.numLoci()-1);
			// decode loci should be fine now
						
			assert (this.decodeLoci.size() > 1);
			this.transitionIndices = new ArrayList<Integer>();
			int maximumStep = 0;
			// now fill the transition indices
			for (int decodeIdx=1; decodeIdx<this.decodeLoci.size(); decodeIdx++) {
				int increment = this.decodeLoci.get(decodeIdx) - this.decodeLoci.get(decodeIdx-1);
				assert ((increment % 2 == 0) || (increment == 1));
				this.transitionIndices.add (increment);
				maximumStep = Math.max (maximumStep, increment);
			}
			// should be fine
			assert (this.transitionIndices.size() == this.decodeLoci.size()-1);
			this.maxStep = maximumStep;
		}


		public static Integer intLogTwo (int integer) {
			assert (integer > 0);
			int number = 1;
			int power = 0;
			while (number <= integer) {
				number *= 2;
				power++;
			}
			return power-1;
		}

		// this parametrization is not nice, but the interfaces don't allow something else
		private boolean isStationaryUniformDistribution (FSAParamSet pSet, int locus, double EPSILON) {
			// the target value
			double uniform = 1d/pSet.numAlleles();
			
			// now see, whether all of them work
			for (int allele=0; allele<pSet.numAlleles(); allele++) {
				if (!isCloseEnough (uniform, pSet.getStationaryProbability(locus, allele), EPSILON)) {
					return false;
				}
			}
			
			// everything fine
			return true;
		}

		// we do this excluding the left one and including the right one
		private Collection<Integer> getPowerOfTwoStepsBetween (int left, int right) {
						
			Collection<Integer> toReturn = new ArrayList<Integer>();
			
			if (left == right) toReturn.add (right);
			
			while (left < right) {
				
				// go right by the biggest power of two that is in the difference
				left += Integer.highestOneBit (right-left);
				
				assert(left <= right);
				
				toReturn.add(left);
			}
			
			return toReturn;
		}

		// this is wrong if theta is changing
		private static boolean checkMutationAndRecombinationRatesAllEqual(FSAParamSet pSet, double epsilon){
			assert (pSet.numLoci() > 1);
			double rho = pSet.getRecombinationRate(0, 1);
			double theta = pSet.getMutationRate(0);
			
			for (int l = 0; l < pSet.numLoci(); l++){
				if (!isCloseEnough(theta, pSet.getMutationRate(l), epsilon)) return false;
			}
			
			for (int l = 0; l < pSet.numLoci() - 1; l++){
				if (!isCloseEnough(rho, pSet.getRecombinationRate(l, l+1), epsilon)) return false;
			}
			
			return true;
		}
		
		private static boolean isCloseEnough(double one, double two, double epsilon){
			return Math.abs(one - two) < epsilon;
		}
		
		
		// this is actually only PIM in very special cases
		// we basically check whether all diagonal elements are the same, and all the off diagonals also
		public static boolean isUniformMutationMatrix (double[][] mutMatrix, double EPSILON){
			
			double diagonal = mutMatrix[0][0];
			double offDiagonal = mutMatrix[0][1];
			
			int numRows = mutMatrix.length;
			
			for (int i = 0; i < numRows; i++){
				if (mutMatrix[i].length != numRows) return false;
				
				for (int j = 0; j < numRows; j++){
					if (i == j && !isCloseEnough(mutMatrix[i][j], diagonal, EPSILON)) return false;
					if (i != j && !isCloseEnough(mutMatrix[i][j], offDiagonal, EPSILON)) return false;
				}
			}
			
			return true;
		}
		
		@Override
		public int numHmmSteps() {
			return this.decodeLoci.size();
		}

		@Override
		public int getLocus(int decodeLocusIdx) {
			return this.decodeLoci.get (decodeLocusIdx);
		}

		@Override
		public int getTransitionIndex (int decodeLocusIdx, boolean reverseDirection) {
			assert (!reverseDirection || (decodeLocusIdx > 0));
			assert (reverseDirection || (decodeLocusIdx < this.transitionIndices.size()));
			return this.transitionIndices.get (decodeLocusIdx - (reverseDirection ? 1 : 0));
		}

		@Override
		public SMCSDTransitionList computeTransitions(UberDemographyCore core, double recoScaling) {
			assert (this.maxStep >= 1);
			// now we have to construct a list of transitions
			// at index k there should actually be the log transitions for a step of length 2^k
			
			// preCompute emissions, which this simple locus skipper does assume that they don't depend on the alleles involved
			double[][][] preEmissionArray = core.getLogEmission(this.mutRate, this.firstEigenMutMatrix);
			assert (core.numDemoStates() == preEmissionArray.length);
			double[] emissions = new double[preEmissionArray.length];
			for (int i=0; i<emissions.length; i++) {
				// TODO maybe assert that they are all the same
				emissions[i] = preEmissionArray[i][0][0];
			}
			
			// new transition list
			SMCSDTransitionList theTransitionList = new SMCSDTransitionList();
			
			// fill it
			
			// first one by hand
			UberDemographyCore.RecoLogProbs recomLogProbs = core.getRecoLogProbs(this.recoRate * recoScaling);	//TODO: will need to do this for each rho encountered
			theTransitionList.setLogNoReco (1, recomLogProbs.logNoReco);
			theTransitionList.setLogReco (1, recomLogProbs.logReco);
			
			// rest via iteration (Chapman-Kolmogoroff)
			for (int k=2; k<=this.maxStep; k *=2) {
				// get the previous ones
				double[] prevNoRecoTransition = theTransitionList.getLogNoReco (k/2);
				double[][] prevRecoTransition = theTransitionList.getLogReco (k/2);
				
				// compute the next ones

				// no reco
				double[] noRecoTransition = new double[prevNoRecoTransition.length];
				for (int intervalDeme=0; intervalDeme<noRecoTransition.length; intervalDeme++) {
					// remmeber, we are in log space, so addition instead of multiplication
					noRecoTransition[intervalDeme] = 2 * prevNoRecoTransition[intervalDeme] + emissions[intervalDeme]; 
				}
				
				// reco
				// [prevIntervalDeme][nextIntervalDeme]
				double[][] recoTransition = new double[prevRecoTransition.length][prevRecoTransition[0].length];
				assert (prevRecoTransition.length == prevRecoTransition[0].length);
				for (int srcIntervalDeme=0; srcIntervalDeme<recoTransition.length; srcIntervalDeme++) {
					for (int dstIntervalDeme=0; dstIntervalDeme<recoTransition[0].length; dstIntervalDeme++) {
						// compute the sum
						LogSum thisSum = new LogSum(recoTransition.length);
						for (int midIntervalDeme=0; midIntervalDeme<recoTransition.length; midIntervalDeme++) {
							// first part
							double firstLogBracket = prevRecoTransition[srcIntervalDeme][midIntervalDeme];
							if (midIntervalDeme == srcIntervalDeme) {
								firstLogBracket = LogSum.computePairLogSum (firstLogBracket, prevNoRecoTransition[srcIntervalDeme]);
							}
							// second part
							double secondLogBracket = prevRecoTransition[midIntervalDeme][dstIntervalDeme];
							if (dstIntervalDeme == midIntervalDeme) {
								secondLogBracket = LogSum.computePairLogSum (secondLogBracket, prevNoRecoTransition[midIntervalDeme]);
							}
							// put together
							double product = firstLogBracket + secondLogBracket + emissions[midIntervalDeme];
							// and put it in the sum
							thisSum.addLogSummand (product);
						}
						
						// correct it, if we have to
						double value = thisSum.retrieveLogSum();
						if (dstIntervalDeme == srcIntervalDeme) {
							// look up the correction
							double correction = noRecoTransition[srcIntervalDeme];
							
							// subtract it
							assert (value >= correction);
							value = LogSum.computePairLogDifference (value, correction);
							
						}
							
						// put value in array
						recoTransition[srcIntervalDeme][dstIntervalDeme] = value;
					}
				}
				
				// and add them
				theTransitionList.setLogNoReco (k, noRecoTransition);
				theTransitionList.setLogReco (k, recoTransition);
			}
			
			// return it
			return theTransitionList;
		}

		@Override
		public double computeLogStationaryProbability (H additionalHaplotype) {
			// in the constructor we made sure that the stationary distribution at each locus is the uniform distribution (comes from the special mutation matrix)
			// thus the stationary distribution is just the per locus stationary (uniform) distribution raised to the appropriate power
			// which is even easier in log-space
			return - this.numLoci * Math.log (this.numAlleles);
		}

		
		@Override
		public double[] getMutationRate(int decodeLocusIdx) {
			return this.mutRate;
		}

		@Override
		public double getLogEmissionProb(double[][] singleBaseEmission, int hmmStep,
				H trunkHap, int trunkPresentDeme, H addHap) {
			int trueLocus = decodeLoci.get(hmmStep);
			return easyGetEmissionProb(singleBaseEmission, trueLocus, trunkHap, addHap);
		}

		@Override
		public void updateEmissionTable(double condProb, double[][] emissionTable,
				int hmmStep, H trunkHap, H addHap, double[][] logEmissionArray) {
			int trueLocus = decodeLoci.get(hmmStep);
			easyUpdateEmissionTable(condProb, emissionTable, trueLocus, trunkHap, addHap);
		}

		
		
		@Override
		public boolean isValidHaplotype (H additionalHaplotype) {
			// is it in the master config?
			return this.validSequences.contains(additionalHaplotype);
		}
		
		@Override
		public boolean isValidConfig (HaplotypeConfiguration<H> config) {
			// go through and see
			for (GeneticTypeMultiplicity<H> hapMult : config) {
				if (!isValidHaplotype(hapMult.geneticType)) return false;
			}
			return true;
		}
		
		protected TreeSet<Integer> getSegregatingSites (Set<H> sequences) {

			TreeSet<Integer> segregatingSites = new TreeSet<Integer>();
						
			// go through
			for (int l=0; l < this.numLoci; l++) {
				
				// it's not segregating if all guys are equal
				Iterator<H> hapIterator = sequences.iterator();
				int firstGuysAllele = hapIterator.next().getAllele(l);
				while (hapIterator.hasNext()) {
					if (firstGuysAllele != hapIterator.next().getAllele(l)) {
						// we found a segregating site
						segregatingSites.add(l);
						break;
					}
				}
			}
			return segregatingSites;
		}
		
		@Override
		public void printInfo (PrintStream outStream) {
			System.out.println ("# Total number loci: " + this.numLoci);
			System.out.println ("# LocsSkipper -- # segregating sites: " + this.numSegSites);
			System.out.println ("# LocsSkipper -- # decode Loci: " + this.decodeLoci.size());
		}
		
		@Override
		public int numNonEmissionLoci() {
			// we don't have explicit emissions for the non-decode loci
			return this.numLoci - this.decodeLoci.size();
		}
		
		private final int numSegSites;
		
		// the decode loci
		final private ArrayList<Integer> decodeLoci;
		
		// and the indices of the steps in between
		final private ArrayList<Integer> transitionIndices;
		final int maxStep;
		
		// the parameters
		private final double recoRate;
		protected double[] mutRate;
		private final MyEigenDecomposition firstEigenMutMatrix;
		
		// we need this for the stationary distribution
		private final int numLoci;
		private final int numAlleles;
		
		private final Set<H> validSequences;

	}
	
	public static class SimpleXLocusSkipper extends SimpleLocusSkipper<HFSAXFullHaplotypeShell> {
		
		public SimpleXLocusSkipper(List<HFSAXFullHaplotypeShell> validSequenceList, EigenParamSet pSet, TreeSet<Integer> additionalDecodeLoci) {
			super(validSequenceList, pSet, additionalDecodeLoci);
		}

		
		
		@Override
		public void updateTheta(
				double[] theta) {
			this.mutRate = Arrays.copyOf(theta, theta.length);
		}



		@Override
		protected TreeSet<Integer> getSegregatingSites(Set<HFSAXFullHaplotypeShell> totalPopulation) {
			return getSNPs(totalPopulation);
		}

		public static TreeSet<Integer> getSNPs(Set<HFSAXFullHaplotypeShell> totalPopulation) {
			FSAReference ref = getReference(totalPopulation);
			assert (ref != null);
			
			TreeSet<Integer> segSites = new TreeSet<Integer>();
			
			for (int i = 0; i < ref.getNumSegLoci(); i++) {
				int toAdd = ref.getSegLocus(i);
				segSites.add(toAdd);
			}
						
			return segSites;
		}

		
		private static FSAReference getReference(Set<HFSAXFullHaplotypeShell> sample) {
			if (sample.size() <= 0) return null;
			
			FSAReference ref = null;
			
			for (HFSAXFullHaplotypeShell H : sample) {
				if (ref == null) {
					ref = H.getReference();
					if (ref == null) return null;
				} else {
					if (ref != H.getReference()) {
						return null;
					}
				}
			}
			
			return ref;
		}	
	}
	
	
	public static class MultiLocusStepHandler<H extends FSAHaplotype> extends HmmStepHandler<H> {

		protected final HashMap<H, Integer> hapToIndexTrunk;
		protected final HashMap<H, Integer> hapToIndexAdditional;
		
		// [addAllele][trunkAllele][addHap][trunkHap][multilocus-block-idx]
		protected final FiveDimShortArray pairAlleleCounts;
		
		private final int numMissingLoci;
		
		protected final int lociPerStep;
		protected final int numSteps;
		
		private final double[] mutRate;
		
		private final int[] hmmstepToRecomRateIdxMap;
		private final List<Double> recomRates;
		
		private final boolean sameSequences;
		
		// the normal constructor
		public MultiLocusStepHandler (List<H> trunkSequences, List<H> additionalSequences, EigenParamSet pSet, int lociPerStep) throws IOException {
			this (trunkSequences, additionalSequences, pSet, lociPerStep, null);
		}
		
		public MultiLocusStepHandler (List<H> trunkSequences, List<H> preAdditionalSequences, EigenParamSet pSet, int lociPerStep, TIntSet partiallyMissing) throws IOException {
			super (pSet);
			
			// are they the same or not
			// only check references for now
			this.sameSequences = (null == preAdditionalSequences);
			
			// deal with the additional sequences
			List<H> additionalSequences = preAdditionalSequences;
			if (this.sameSequences) additionalSequences = trunkSequences;
			
			// get some indices (fill trunk first)
			this.hapToIndexTrunk = new HashMap<H, Integer>();
			for (int i = 0; i < trunkSequences.size(); i++) {
				hapToIndexTrunk.put(trunkSequences.get(i), i);
			}
			// now see about the additional sequences
			if (this.sameSequences) {
				this.hapToIndexAdditional = this.hapToIndexTrunk;
			}
			else {
				// make one for the additional ones
				this.hapToIndexAdditional = new HashMap<H, Integer>();
				for (int i = 0; i < additionalSequences.size(); i++) {
					hapToIndexAdditional.put(additionalSequences.get(i), i);
				}

			}
			// how many site should be lumped together
			this.lociPerStep = lociPerStep;
			
			// how many multiloci do we have
			this.numSteps = pSet.numLoci() / lociPerStep + (pSet.numLoci() % lociPerStep == 0 ? 0 : 1);
			
			
			// for some stuff we might need a combined list
			List<H> combinedList = trunkSequences;
			if (!this.sameSequences) {
				combinedList = new ArrayList<H>(trunkSequences);
				combinedList.addAll(additionalSequences);
			}
			
			
			// get some preSegregating
			TIntSet preSegSites = null;
			TIntSet segSites = null;
			// container [allele][multiLocus] (everything should be intialized to zero)
			short[][] nonSegregatingAlleles = null;
			// do we have a reference ? (for now this is the decisive thing)
			if (combinedList.get(0) instanceof FSAXHaplotype) {
				// we have a reference, so get the segregating sites
				FSAReference reference = ((FSAXHaplotype)combinedList.get(0)).getReference();
				// in a set
				preSegSites = new TIntHashSet ();
				for (int s : reference.getSegSites()) {
					preSegSites.add(s);
				}

				// get some things going for later
				// also set up the container
				nonSegregatingAlleles = new short[this.pSet.numAlleles()][this.numSteps];
				segSites = new TIntHashSet();
			}

//			long startTime = System.currentTimeMillis();
			// figure out how many non-missing bases
			// this handler: one missing in any haplotype immediately declares the whole locus as missing
			// (although we have to account for the other handler)
//			TIntSet missingLoci = new TIntHashSet ();
			// can only take int as size
			BitSet missingLoci = new BitSet (this.pSet.numLoci());

			// loop through
			for (int locus = 0; locus < this.pSet.numLoci(); locus++) {
				
				// everybody starts out not missing
				missingLoci.set (locus, false);

				// if not segregating, we only need to look at the first guy (reference)
				if (preSegSites != null && !preSegSites.contains(locus)) {
					if (combinedList.get(0).getAllele(locus) == ReadSequences.MISSING_ALLELE) {
						missingLoci.set (locus, true);
					}
					// next one
					continue;
				}
				
				// now it is in preSegrating
				
				// count how many missing at this locus
				int numMissing = 0;
				for (H sequence : combinedList) {
					// how many missing at this locus
					if (sequence.getAllele(locus) == ReadSequences.MISSING_ALLELE) {
						numMissing += 1;
					}
				}
				
				// everthing is fine if no one is missing
				if (numMissing == 0) {
					// not missing, only segregating
					if (segSites != null) segSites.add(locus);
				}
				// what if there is something missing
				else{
					if ((numMissing != combinedList.size()) && (partiallyMissing != null)) {
						// yes, we care for partially
						partiallyMissing.add(locus);
						// don't add segregating partially missing
					}
					else {
						// just count as missing
						missingLoci.set (locus, true);
						// this does not even matter
						if (segSites != null) segSites.add(locus);
					}
				}
			}
			// set the non-missing bases
			// works with set
//			this.numMissingLoci = missingLoci.size();
			// need to change for bitset
			this.numMissingLoci = missingLoci.cardinality();

			// just for convenience we also want the smaller resolution
			// find the differences
			//[addAllele][trunkAllele][addHap][trunkHap][multilocus-block-idx]
			try {
				this.pairAlleleCounts = new FiveDimShortArray (this.pSet.numAlleles(), this.pSet.numAlleles(), additionalSequences.size(), trunkSequences.size(), numSteps);
			}
			catch (IOException e) {
				throw new IOException (e.getMessage() + " Try using shorter sequences, less individuals, a larger binSize, or fewer alleles.");
			}

			// go through
			for (int locus = 0; locus < this.pSet.numLoci(); locus++) {
				// where we at (the floor)
				int currStep = locus / lociPerStep;
				
				// don't use it if it is missing (don't have to)
				// partially missing will be treated elsewhere too
				if (missingLoci.get(locus) || ((partiallyMissing != null) && (partiallyMissing.contains(locus)))) {
					continue;
				}
				
				// if it is non-segregating, we can do this faster too
				if ((segSites != null) && !segSites.contains(locus)) {
					// remember the allele
					nonSegregatingAlleles[trunkSequences.get(0).getAllele(locus)][currStep] += 1;
					// important, go out of the loop
					continue;
				}
				
				// get the allele counts between the different haplotypes
				// this does too much, but ok for now
				for (int addIdx = 0; addIdx < additionalSequences.size(); addIdx++) {
					
					Integer trunkStart = addIdx;
					if (!this.sameSequences) trunkStart = 0;
					
					for (int trunkIdx = trunkStart; trunkIdx < trunkSequences.size(); trunkIdx++) {
						// we need to do do self-comparing
						// also, see what case we are in
						int addAllele = additionalSequences.get(addIdx).getAllele(locus);
						int trunkAllele = trunkSequences.get(trunkIdx).getAllele(locus);
						// nobody missing
						//[addAllele][trunkAllele][addHap][trunkHap][multilocus-block-idx]
						this.pairAlleleCounts.adjust (addAllele, trunkAllele, addIdx, trunkIdx, currStep, (short) 1);
						// only do the mirroring if same sequences
						// otherwise we should already iterate through everything
						if (this.sameSequences && (addIdx != trunkIdx)) this.pairAlleleCounts.adjust (trunkAllele, addAllele, trunkIdx, addIdx, currStep, (short) 1);
					}
				}
			}
			
			// if we have info about non-segregating sites, we should add it
			// but only in the allele pair count, i think
			if (nonSegregatingAlleles != null) {
				// for each pair of haplotypes
				for (int addIdx = 0; addIdx < additionalSequences.size(); addIdx++) {
					for (int trunkIdx = 0; trunkIdx < trunkSequences.size(); trunkIdx++) {
						// we need to do self comparing
						// this has to be done for each lump
						for (int allele=0; allele<this.pSet.numAlleles(); allele++) {
							for (int currStep=0; currStep<this.numSteps; currStep++) {
								this.pairAlleleCounts.adjust (allele, allele, addIdx, trunkIdx, currStep, nonSegregatingAlleles[allele][currStep]);
							}
						}
					}
				}
			}

			// only one reco and one mut rate
			// with the new type for all parameters same everywhere, this is given
//			assert (SimpleLocusSkipper.checkMutationAndRecombinationRatesAllEqual(pSet, UberDemographyCore.ONELOCUS_EPSILON));

			// there can be only one
			mutRate= new double[1];
			mutRate[0] = this.pSet.getMutationRate(0);
			
			// same here, store special, for now just multiply
			recomRates = new ArrayList<Double>();
			recomRates.add(this.pSet.getRecombinationRate(0, 1) * this.lociPerStep);
			
			// fake the map
			this.hmmstepToRecomRateIdxMap = new int[this.numSteps];
			for (int l=0; l<this.hmmstepToRecomRateIdxMap.length; l++) {
				this.hmmstepToRecomRateIdxMap[l] = 0;
			}
			
			// everything should be done now
//			long elapsedTime = System.currentTimeMillis() - startTime;
//			MetaOptimization.synchronizedPrintln ("missing time: " + elapsedTime);
		}
		
		@Override
		public double computeLogStationaryProbability (H haplotype) {
			// did we already compute it?
			Double value = this.stationaryMap.get(haplotype);
			
			if (value == null) {
				// apparently not, so compute it

				int hapIdx = this.hapToIndexTrunk.get (haplotype);

				// here we have to assume that all mutations are equal along the sequence
				// with the new type for all parameters same everywhere, this is given
//				assert (SMCSDemo.checkMutationMatricesEqual (pSet, UberDemographyCore.ONELOCUS_EPSILON, UberDemographyCore.RELATIVE_ERROR));
				double[] logStationaryProb = new double[pSet.numAlleles()];
				for (int a=0; a<logStationaryProb.length; a++) {
					logStationaryProb[a] = Math.log (this.pSet.getStationaryProbability(0, a));
				}
				
				// but everything in log scale
				value = 0d;
				// go though meta-loci, last size is the number of loci (due to stupid java)
				for (int l=0; l<this.pairAlleleCounts.getFifthSize(); l++) {
					// go through the alleles, just ignoring missing handles them correctly
					for (int a=0; a<pSet.numAlleles(); a++) {
						// how many a's?
						//[addAllele][trunkAllele][addHap][trunkHap][multilocus-block-idx]
						int numA = this.pairAlleleCounts.get(a, a, hapIdx, hapIdx, l);
						value += logStationaryProb[a] * numA;
					}
				}
				// should be it
				
				// put it in the map
				stationaryMap.put(haplotype, value);
			}
			
			// return it
			return value;
		}

		
		@Override
		public int numHmmSteps() {
			return this.numSteps;
		}

		@Override
		public int numNonEmissionLoci() {
			// doesn't emit at them missing loci
			return this.numMissingLoci;
		}
		
		@Override
		public int getTransitionIndex(int decodeLocusIdx, boolean reverseDirection) {
			assert (!reverseDirection || (decodeLocusIdx > 0));
			assert (reverseDirection || (decodeLocusIdx < this.hmmstepToRecomRateIdxMap.length));
			return this.hmmstepToRecomRateIdxMap[decodeLocusIdx - (reverseDirection ? 1 : 0)];
		}

		@Override
		public double[] getMutationRate(int decodeLocusIdx) {
			return mutRate;
		}

		@Override
		public SMCSDTransitionList computeTransitions(UberDemographyCore core, double recoScaling) {
			return easyComputeTransitions(core, recomRates, recoScaling);
		}
		
		@Override
		public double getLogEmissionProb(double[][] singleBaseLogEmission, int hmmStep, H trunkHap, int trunkPresentDeme, H addHap) {
			Double toReturn = computeEmission (null, null, singleBaseLogEmission, hmmStep, trunkHap, addHap);
			assert (toReturn != null);
			return toReturn;
		}

		@Override
		public void updateEmissionTable(double condProb, double[][] emissionTable, int hmmStep, H trunkHap, H addHap, double[][] logEmissionArray) {
			Double toReturn = computeEmission (condProb, emissionTable, logEmissionArray, hmmStep, trunkHap, addHap);
			assert (toReturn == null);
		}

		private Double computeEmission (Double condProb, double[][] emissionTable, double[][] singleBaseLogEmission, int hmmStep, H trunkHap, H addHap) {
			// this should never happen, otherwise it's used incorrectly
			assert ((condProb == null) == (emissionTable == null));
			
			// now do your thing
			int trunkIdx = this.hapToIndexTrunk.get(trunkHap);
			int addIdx = this.hapToIndexAdditional.get(addHap);
			
			int numAlleles = singleBaseLogEmission.length;

			// just have it
			double toReturn = 0d;
			
			// loop over actual alleles
			for (int addAllele = 0; addAllele < numAlleles; addAllele++) {
				for (int trunkAllele = 0; trunkAllele < numAlleles; trunkAllele++) {
					// either just compute the emission
					if (emissionTable == null) {
						// this if-check properly handles case where singleBaseLogEmission = infinity
						//[addAllele][trunkAllele][addHap][trunkHap][multilocus-block-idx]
						if (this.pairAlleleCounts.get (addAllele, trunkAllele, addIdx, trunkIdx, hmmStep) != 0) {
							toReturn += singleBaseLogEmission[trunkAllele][addAllele] * this.pairAlleleCounts.get (addAllele, trunkAllele, addIdx, trunkIdx, hmmStep);
						}
					}
					// or update the emission table
					else {
						//[addAllele][trunkAllele][addHap][trunkHap][multilocus-block-idx]
						emissionTable[trunkAllele][addAllele] += condProb * this.pairAlleleCounts.get (addAllele, trunkAllele, addIdx, trunkIdx, hmmStep);
					}
				}
			}
			
			// only return something if we wanted to compute a single probability
			if (emissionTable == null) {
				return toReturn;
			}
			else {
				return null;
			}
		}
		
		@Override
		public boolean isValidHaplotype(H additionalHaplotype) {
			if (this.sameSequences) {
				return hapToIndexAdditional.containsKey(additionalHaplotype);
			}
			else {
				return (hapToIndexAdditional.containsKey(additionalHaplotype) || hapToIndexTrunk.containsKey(additionalHaplotype));
			}
		}

		@Override
		public boolean isValidConfig(HaplotypeConfiguration<H> config) {
			// go through and see
			for (GeneticTypeMultiplicity<H> hapMult : config) {
				if (!isValidHaplotype(hapMult.geneticType)) return false;
			}
			return true;
		}

		@Override
		public void printInfo(PrintStream outStream) {
			System.out.println("# Multilocus-block HMM type");
			System.out.println("# Num multilocus-blocks: " + this.numSteps);
			System.out.println("# Num loci per block: " + this.lociPerStep);
		}
		
	}

	
	public static class MultiLocusCleverImputeStepHandler<H extends FSAHaplotype> extends MultiLocusStepHandler<H> {

		
		private final HmmStepHandler.SingleStepMissingAllelesTransitionMap.CoreCache<H> coreCache;
		private final DemoConfiguration<H> configs;

		public MultiLocusCleverImputeStepHandler(H addHap,
				EigenParamSet pSet, int lociPerStep, HmmStepHandler.SingleStepMissingAllelesTransitionMap.CoreCache<H> coreCache,
				DemoConfiguration<H> configs) throws IOException {
			super(getSequenceList(configs, addHap), null, pSet, lociPerStep);
			
			this.coreCache = coreCache;
			this.configs = configs;
		}

		@Override
		public double getLogEmissionProb(double[][] singleBaseLogEmission,
				int hmmStep, H trunkHap, int trunkPresentDeme, H addHap) {
			double toReturn = 0d;
			for (int locus = hmmStep * this.lociPerStep; locus < Math.min((hmmStep+1) * lociPerStep, this.pSet.numLoci()); locus++) {
				double[] imputedAlleles = SingleStepMissingAllelesTransitionMap.getTrunkHapAlleleProbs(locus, trunkHap, trunkPresentDeme, configs, coreCache, pSet);
				
				toReturn += SingleStepMissingAllelesTransitionMap.oneLocusEmissionProb(singleBaseLogEmission, locus, imputedAlleles, trunkPresentDeme, addHap);
			}
			return toReturn;
		}
		
		public static <H extends FSAHaplotype> List<H> getSequenceList(DemoConfiguration<H> configs, H addHap) {
			HaplotypeConfiguration<H> totalConfig = configs.getTotalPopulation();
			
			List<H> toReturn = new ArrayList<H>();
			for (GeneticTypeMultiplicity<H> hapMult : totalConfig) {
				if (hapMult.geneticType == addHap) continue;
				toReturn.add(hapMult.geneticType);
			}
			toReturn.add(addHap);
			
			
			return toReturn;
		}
		
	}
	
	
	// now some inheritance
	public static class MultiLocusStepHandlerSuboptimalMissing<H extends FSAHaplotype> extends MultiLocusStepHandler<H> implements HasPartiallyMissing{

		// TODO REORDER THESE TOO
		// we only need the partially missing stuff
		private final TIntSet partiallyMissing;
		//[multilocus-block-idx][hap1][hap2][allele]
		private final int[][][][] firstMissing;
		//[multilocus-block-idx][hap1][hap2]
		private final int[][][] bothMissing;
		
		// stationary is needed here
		private final double[] alleleLogStationaryDistn;

		private final int numPartiallyMissing;

		// object oriented woodoo
		public MultiLocusStepHandlerSuboptimalMissing (List<H> trunkSequences, List<H> additionalSequences, EigenParamSet pSet, int lociPerStep) throws IOException {
			this (trunkSequences, additionalSequences, pSet, lociPerStep, new TIntHashSet());
		}

		private MultiLocusStepHandlerSuboptimalMissing (List<H> trunkSequences, List<H> additionalSequences, EigenParamSet pSet, int lociPerStep, TIntSet partiallyMissing) throws IOException {
			super (trunkSequences, additionalSequences, pSet, lociPerStep, partiallyMissing);
			
			
			
			
			
			
			// TODO REORDER THE ARRAYS to have the largest number last
			// also we have to see about the stuff with distinguishing between trunk sequences and additional sequences
			if (true) throw new IOException("handling missing alleles with the stationary distribution currently not implemented for multi locus stephandler.");
			
			// something, but I'm sure it won't work properly
			List<H> allSequences = new ArrayList<H>(trunkSequences);
			allSequences.addAll(additionalSequences);
			
			
			
			
			
			this.partiallyMissing = partiallyMissing;

			this.numPartiallyMissing = partiallyMissing.size();
			
			// everything that does not involve partially missing sites should have been done in the super constructor
						
			// the stationary, first locus should be fine
			double[] stat = StationaryDistribution.getStationaryDistribution(this.pSet.getMutationMatrix(0).getArray());
			this.alleleLogStationaryDistn = new double[stat.length];
			for (int i=0; i<stat.length; i++) {
				this.alleleLogStationaryDistn[i] = Math.log(stat[i]);
			}

			// we now want to deal with partially missing in a different way
			// partially missing should have been filled in the super constructor
			// find the differences (that involve partially missing)
			this.firstMissing = new int[this.numSteps][allSequences.size()][allSequences.size()][this.pSet.numAlleles()];
			this.bothMissing = new int[this.numSteps][allSequences.size()][allSequences.size()];
			
			
			// go through
			for (int locus : this.partiallyMissing.toArray()) {
				// where we at (the floor)
				int currStep = locus / lociPerStep;
								
				// get the allele counts between the different haplotypes
				for (int i = 0; i < allSequences.size(); i++) {
					for (int j = 0; j < allSequences.size(); j++) {
						// also, see what case we are in
						int hapOneAllele = allSequences.get(i).getAllele(locus);
						int hapTwoAllele = allSequences.get(j).getAllele(locus);
						if (hapOneAllele == ReadSequences.MISSING_ALLELE) {
							// one is missing
							if (hapTwoAllele == ReadSequences.MISSING_ALLELE) {
								// second also
								this.bothMissing[currStep][i][j] += 1;
							}
							else {
								// only first missing
								this.firstMissing[currStep][i][j][hapTwoAllele] += 1;
							}
						}
						else {
							// we actually don't need to care about allele two missing, cause that's taken care of at another point in the loop
							if (hapTwoAllele != ReadSequences.MISSING_ALLELE) {
								// oh, in this case we actually do have to invoke the pairAlleleCounts
								this.pairAlleleCounts.adjust (hapOneAllele, hapTwoAllele, i, j, currStep, (short) 1);
							}
						}
					}
				}
			}
			// everything should be done now
		}
		
		

		
		@Override
		public double getLogEmissionProb(double[][] singleBaseLogEmission, int hmmStep, H trunkHap, int trunkPresentDeme, H addHap) {
			// call the super thing
			double superResult = super.getLogEmissionProb(singleBaseLogEmission, hmmStep, trunkHap, trunkPresentDeme, addHap);
			// then call our method
			Double thisResult = computeEmission (null, null, singleBaseLogEmission, hmmStep, trunkHap, addHap);
			assert (thisResult != null);
			// return them both
			return (superResult + thisResult);
		}

		@Override
		public void updateEmissionTable(double condProb, double[][] emissionTable, int hmmStep, H trunkHap, H addHap, double[][] logEmissionArray) {
			// let the super method do it's thing
			super.updateEmissionTable(condProb, emissionTable, hmmStep, trunkHap, addHap, logEmissionArray);
			// then we need to use this method to update the emission table further
			Double toReturn = computeEmission (condProb, emissionTable, logEmissionArray, hmmStep, trunkHap, addHap);
			// and we should be fine
			assert (toReturn == null);
		}

		private Double computeEmission (Double condProb, double[][] emissionTable, double[][] singleBaseLogEmission, int hmmStep, H trunkHap, H addHap) {
			// this should never happen, otherwise it's used incorrectly
			assert ((condProb == null) == (emissionTable == null));
			
			// now do your thing
			int trunkIdx = this.hapToIndexTrunk.get(trunkHap);
			int addIdx = this.hapToIndexAdditional.get(addHap);
			
			
			assert false;
			// TODO some of the indices might have to be changed
			// the whole class doesn't work properly right now
			
			
			int numAlleles = singleBaseLogEmission.length;
			assert (numAlleles == this.alleleLogStationaryDistn.length);

			// just have it
			double toReturn = 0d;
			
			// loop over actual alleles
			// but only do the things for the partially missing
			for (int addAllele = 0; addAllele < numAlleles; addAllele++) {
				
				// if we emit this allele from nothing, we have to multiply by the stationary distribution
				// if we emit nothing, we don't have to modify the log prob at all
				if (emissionTable == null) {
					// don't do nothing if this wasn't observed
					if (this.firstMissing[hmmStep][trunkIdx][addIdx][addAllele] != 0) {
						toReturn += this.alleleLogStationaryDistn[addAllele] * this.firstMissing[hmmStep][trunkIdx][addIdx][addAllele];
					}
				}
				// if we emit nothing, no changes to emission table
				// if we emit this allele from nothing, we have to update the emission table based on the posterior probability
				else {
					// don't do nothing if we need not
					if (this.firstMissing[hmmStep][trunkIdx][addIdx][addAllele] == 0) continue;
					// compute posterior
					double[] post = getStupidStationaryTrunkAllelePosterior (addAllele, condProb, this.alleleLogStationaryDistn, singleBaseLogEmission);
					
					for (int k=0; k<post.length; k++) {
						emissionTable[k][addAllele] += post[k] * this.firstMissing[hmmStep][trunkIdx][addIdx][addAllele];
					}
				}
			}
			
			// only return something if we wanted to compute a single probability
			if (emissionTable == null) {
				return toReturn;
			}
			else {
				return null;
			}
		}

		@Override
		// we need this here
		public int numPartiallyMissing() {
			return this.numPartiallyMissing;
		}

	}

	
	public static double[] getStupidStationaryTrunkAllelePosterior (int addAllele, double condProb, double[] logStationary, double[][] logEmissionProbs) {
		
		int numAlleles = logStationary.length;
		assert ((0 <= addAllele) && (addAllele < numAlleles));
		
		// conditional probability of hidden state and trunk allele, given the data
		double[] alleleCondProbs = new double[numAlleles];
		
		double sum = 0d;
		for (int trunkAllele = 0; trunkAllele < numAlleles; trunkAllele++) {
			alleleCondProbs[trunkAllele] = Math.exp (logStationary[trunkAllele] + logEmissionProbs[trunkAllele][addAllele]);
			sum += alleleCondProbs[trunkAllele];
		}
		
		assert Math.abs(sum - Math.exp (logStationary[addAllele])) < UberDemographyCore.ONELOCUS_EPSILON;
		
		// alleleCondProbs should sum to condProb, but it doesnt because of constant of proportionality
		// so correct
		for (int trunkAllele = 0; trunkAllele < numAlleles; trunkAllele++) {
			alleleCondProbs[trunkAllele] *= condProb / sum;
		}

		// give them away now
		return alleleCondProbs;
	}

}
