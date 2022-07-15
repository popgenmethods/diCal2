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

import java.util.Iterator;
import java.util.LinkedHashMap;

import edu.berkeley.diCal2.haplotype.GeneticType;
import edu.berkeley.diCal2.haplotype.GeneticType.GeneticTypeMultiplicity;

// NOTE WELL: template parameters of this class are used in a HashMap, and so must either implement hashCode()
//  and equals() or otherwise ensure that hashing works properly (e.g. private constructor + factory)

public class HaplotypeConfiguration<H extends GeneticType> implements Iterable<GeneticTypeMultiplicity<H>> {
	// the configuration
	private final LinkedHashMap<H, GeneticTypeMultiplicity<H>> typeConfiguration;
	
	// cached information about the configuration
	private int numGeneticTypes;
	private int distinctGeneticTypes;

	// what can we expect
	public final int numLoci;
	public final int numAlleles;
	
	public HaplotypeConfiguration (int numLoci, int numAlleles) {
		this.typeConfiguration = new LinkedHashMap<H, GeneticTypeMultiplicity<H>>();
		
		this.numGeneticTypes = 0;
		this.distinctGeneticTypes = 0;
		
		this.numLoci = numLoci;
		this.numAlleles = numAlleles;
	}

	public int numLoci() { return numLoci; }
	public int numAlleles() { return numAlleles; }
	
	public final int totalGeneticTypes() {
		return numGeneticTypes;
	}
	
	public final int distinctGeneticTypes() {
		return distinctGeneticTypes;
	}

	public void adjustType(H type, int adjustment) {
		GeneticTypeMultiplicity<H> geneMultiplicity = typeConfiguration.get(type);
		
		if (geneMultiplicity == null) {
			if (adjustment != 0) {
				typeConfiguration.put(type, new GeneticTypeMultiplicity<H>(type, adjustment));
				distinctGeneticTypes++;
			}
		} else {
			int newMultiplicity = geneMultiplicity.multiplicity + adjustment;
			if (newMultiplicity != 0) {
				geneMultiplicity.multiplicity = newMultiplicity;
			} else {
				typeConfiguration.remove(type);
				distinctGeneticTypes--;
			}
		}
		
		numGeneticTypes += adjustment;
		assert (numGeneticTypes >= 0);
	}

	public final void addTypes(Iterable<GeneticTypeMultiplicity<H>> configuration) {
		for (GeneticTypeMultiplicity<H> typeMultiplicity : configuration) {
			adjustType(typeMultiplicity.geneticType, typeMultiplicity.multiplicity);
		}
	}

	public final int getMultiplicity(H type) {
		GeneticTypeMultiplicity<H> typeMultiplicity = typeConfiguration.get(type);
		return typeMultiplicity == null ? 0 : typeMultiplicity.multiplicity;
	}
	
	@Override
	public final Iterator<GeneticTypeMultiplicity<H>> iterator() {
		return typeConfiguration.values().iterator();
	}

	@Override
	public String toString() {
		StringBuffer sb = new StringBuffer();
		
		for (GeneticTypeMultiplicity<H> typeMultiplicity : this) {
			sb.append(typeMultiplicity.geneticType + ":" + typeMultiplicity.multiplicity + "\n");
		}
		
		return sb.toString();
	}
}