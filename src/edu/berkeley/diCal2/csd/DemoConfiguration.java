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

import java.io.BufferedReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import edu.berkeley.diCal2.haplotype.GeneticType;
import edu.berkeley.diCal2.haplotype.GeneticType.GeneticTypeMultiplicity;
import edu.berkeley.diCal2.maximum_likelihood.DiCalParamSet.DemeHapPair;

public class DemoConfiguration<H extends GeneticType> {
	/// list of configurations (not so nice, cause configurations might have different allele numbers)
	private ArrayList<HaplotypeConfiguration<H>> populations;
	private final int numLoci;
	private final int numAlleles;
	private final int numDemes;
	
	
	/// construct one with a given number of empty demes
	public DemoConfiguration (int numDemes, int numLoci, int numAlleles) {
		this(numDemes, numLoci, numAlleles, new ArrayList<H>(), new int[0][]);
	}
	
	// constructs a filled DemoConfiguration from the haplotypes and config Info
	public DemoConfiguration (ConfigInfo configInfo, List<H> haplotypes) {
		this(configInfo.numDemes, configInfo.numLoci, configInfo.numAlleles, haplotypes, configInfo.multiplicities.toArray(new int[][] {}));
	}
	
	// constructs a filled DemoConfiguration
	public DemoConfiguration (int numDemes, int numLoci, int numAlleles, List<H> haplotypes, int[][] multiplicities) {
		
		// save some things
		this.numLoci = numLoci;
		this.numAlleles = numAlleles;
		this.numDemes = numDemes;
		
		// make it empty (might not be necessary)
		populations = new ArrayList<HaplotypeConfiguration<H>> ();
		
		// fill it with enough empty ones
		for (int i=0; i<numDemes; i++) {
			populations.add (new HaplotypeConfiguration<H>(numLoci, numAlleles));
		}
		
		assert (haplotypes.size() == multiplicities.length);
		
		for (int hapIdx = 0; hapIdx < haplotypes.size(); hapIdx++) {
			H currHap = haplotypes.get(hapIdx);
			int[] currMult = multiplicities[hapIdx];
			assert (currMult.length == this.numDemes);
			
			for (int deme = 0; deme < this.numDemes; deme++) {
				this.adjustHaplotypeMultiplicity(deme, currHap, currMult[deme]);
			}
		}
	}
	
	
	
	public DemoConfiguration(int numDemes, int numLoci, int numAlleles, List<DemeHapPair<H>> orderedSampleConfig) {
		// save some stuff
		this.numLoci = numLoci;
		this.numAlleles = numAlleles;
		this.numDemes = numDemes;
		
		// make it empty (might not be necessary)
		populations = new ArrayList<HaplotypeConfiguration<H>> ();
		
		// fill it with enough empty guys
		for (int i=0; i<numDemes; i++) {
			populations.add (new HaplotypeConfiguration<H>(numLoci, numAlleles));
		}
		
		for (DemeHapPair<H> demeHap : orderedSampleConfig) {
			this.adjustHaplotypeMultiplicity(demeHap.demeNumber, demeHap.haplotype, 1);
		}
	}

	// need to make a copy 
	public static <H extends GeneticType> DemoConfiguration<H> copyConfiguration (DemoConfiguration<H> toCopy) {
		
		// do some of the copying
		DemoConfiguration<H> returnConfig = new DemoConfiguration<H> (toCopy.getNumberDemes(), toCopy.numLoci(), toCopy.numAlleles());
		
		// just copy every population
		for (int i=0; i<toCopy.getNumberDemes(); i++) {
			returnConfig.addHaplotypes(i, toCopy.getPopulation(i));
		}
		
		// return it
		return returnConfig;
	}
	
	public int numLoci () {
		for (HaplotypeConfiguration<H> pop : this.populations) {
			assert (this.numLoci == pop.numLoci());
		}
		return this.numLoci;
	}
	
	public int numAlleles () {
		for (HaplotypeConfiguration<H> pop : this.populations) {
			assert (this.numAlleles == pop.numAlleles());
		}
		return numAlleles;
	}
	
	public boolean isGhost (int demeNumber) {
		return getPopulation(demeNumber).totalGeneticTypes() == 0;
	}
	
	/// how many demes?
	public int getNumberDemes () {
		assert (this.numDemes == populations.size());
		return populations.size();
	}
	
	/// get a population by deme
	public HaplotypeConfiguration<H> getPopulation (int demeNumber) {
		// check number
		assert (0 <= demeNumber && demeNumber < getNumberDemes());
		// give it if ok
		return populations.get (demeNumber);
	}
	
	/// set a population in a given deme
	public void setPopulation (int demeNumber, HaplotypeConfiguration<H> population) {
		// check number
		assert (0 <= demeNumber && demeNumber < getNumberDemes());
		// set it, if ok
		populations.set (demeNumber, population);
	}
	
	/// adjust a haplotype multiplicity in a given population
	public void adjustHaplotypeMultiplicity (int demeNumber, H haplotype, int adjustment) {
		// check number
		assert (0 <= demeNumber && demeNumber < getNumberDemes());
		// and adjust
		(populations.get (demeNumber)).adjustType (haplotype, adjustment);
	}

	/// add some haplotypes to a given population
	public void addHaplotypes (int demeNumber, Iterable<GeneticTypeMultiplicity<H>> configuration) {
		// check number
		assert (0 <= demeNumber && demeNumber < getNumberDemes());
		// and adjust
		(populations.get (demeNumber)).addTypes (configuration);
	}
	
	/// give some nice output
	@Override
	public String toString () {
		String result = "";
		// go through guys
		for (int i=0; i<getNumberDemes(); i++) {
			// say what deme we are in
			result += "==== " + i +  " ====\n";
			result += populations.get(i).toString();
		}
		return result;
	}

	public HaplotypeConfiguration<H> getTotalPopulation () {
		// get an index set with all demes
		int[] indexSet = new int[this.getNumberDemes()];
		for (int i=0; i<indexSet.length; i++) {
			indexSet[i] = i;
		}
		// and return the combined subpopulations
		return this.getCombinedSubpopulation(indexSet);
	}
	
	public HaplotypeConfiguration<H> getCombinedSubpopulation(int[] subPops) {
		// put everything in this one
		HaplotypeConfiguration<H> thePop = new HaplotypeConfiguration<H> (this.populations.get(0).numLoci(), this.populations.get(0).numAlleles());
		
		// go through given numbers (doubles are ok)
		for (int gamma : subPops) {
			// make sure its actually a real deme
			assert ((gamma >= 0) && (gamma < this.getNumberDemes()));
			
			// add the multiplicities of the requested populations
			for (GeneticTypeMultiplicity<H> hapMult : this.populations.get(gamma)) {
				thePop.adjustType (hapMult.geneticType, hapMult.multiplicity);
			}
		}
		
		// give it away now
		return thePop;
	}
	
	public int[] getSampleSizes() {
		int[] n= new int[getNumberDemes()];
		for (int i=0; i < populations.size() ; i++) {
			n[i]= populations.get(i).totalGeneticTypes();
		}
		return n;
	}
	
	
	public int getTotalGeneticTypes() {
		int n=0;
		for (HaplotypeConfiguration<H> thisPop : populations) {
			n += thisPop.totalGeneticTypes();
		}
		return n;
	}

	public int getDistinctGeneticTypes() {
		int n=0;
		for (HaplotypeConfiguration<H> thisPop : populations) {
			n += thisPop.distinctGeneticTypes();
		}
		return n;
	}
	
	public static class ConfigInfo {
		public List<int[]> multiplicities;
		public int numDemes;
		public Integer numLoci;
		public int numAlleles;
		
		public ConfigInfo(List<int[]> multiplicities, int numDemes, Integer numLoci, int numAlleles) {
			this.multiplicities = multiplicities;
			this.numDemes = numDemes;
			this.numLoci = numLoci;
			this.numAlleles = numAlleles;
		}
	}
	
	public static ConfigInfo readConfigInfo(BufferedReader bufferedConfigReader) throws IOException {
		
		// go through lines
		int numDemes=0, numAlleles=0;
		long numLoci=0;
		ArrayList<int[]> multiplicites = new ArrayList<int[]>();
		String line = null;
		boolean expectHeader = true;
		while ((line = bufferedConfigReader.readLine()) != null) {
			// always ignore comment and empty lines
			if (line.startsWith("#") || line.trim().equals("")) continue;
			
			// maybe header?
			if (expectHeader) {
				// first line in numSites, numAlleles and numDemes
				String[] fields = line.split("(\t| )+");
				if (fields.length != 3) throw new IOException("First line (header) of config file not formatted correctly. Expected: <numLoci> <numAlleles> <numDemes>");
				
				// get data
				numLoci = Long.parseLong(fields[0]);
				numAlleles = Integer.parseInt(fields[1]);
				numDemes = Integer.parseInt(fields[2]);
				
				// no header anymore
				expectHeader = false;
			}
			else {
				// body
				// which should be some multiplicities for a haplotype
				String[] fields = line.split("(\t| )+");
				if (fields.length != numDemes) throw new IOException("Multiplicity line in config file has invalid number of columns.");
				
				// convert it
				int[] mult = new int[numDemes];
				for (int i=0; i<fields.length; i++) {
					mult[i] = Integer.parseInt(fields[i]);
				}
				
				// and add it
				multiplicites.add (mult);
			}
		}
		// should be fine for now

		// see about whether numLoci is too large
		// TODO: real solution is to make numLoci a long in config, but there is too many places that use int
		Integer realNumLoci = null;
		if (numLoci < Integer.MAX_VALUE) {
			realNumLoci = new Integer ((int) numLoci);
		}
		
		// give it away now
		return new ConfigInfo(multiplicites, numDemes, realNumLoci, numAlleles);
	}
	
	// removes all haplotypes of multiplicity 0 from the lists
	public static <H extends GeneticType> void removeMissingHaplotypes (ConfigInfo configInfo, List<H> hapList) {
		assert (configInfo.multiplicities.size() == hapList.size());
		
		for (int idx = hapList.size() - 1; idx >= 0; idx--) {
			boolean allZero = true;
			for (int mult : configInfo.multiplicities.get(idx)) allZero &= (mult == 0);
			if (allZero) {
				hapList.remove(idx);
				configInfo.multiplicities.remove(idx);
			}
		}
	}
	}
