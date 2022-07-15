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


public interface FSAXTypeConfig<H extends FSAXHaplotype> extends GeneticTypeConfig<H> {
	public int numLoci();
	public int numAlleles();
	public FSAReference getReference();
	
	public interface FSAReference {
		public boolean locusSegregating(int locus);
		public int getNonSegAllele(int locus);
		
		public int getLocalIdx(int segLocus);
		public int getSegLocus(int localIdx);
		
		public int getNumLoci();
		public int getNumSegLoci();
		
		public int[] getSegSites();
	}
	
	public interface FSAXConfigConstructor<H extends FSAXHaplotype, C extends FSAXTypeConfig<H>> {
		public C newFSAXConfig(int numLoci, int numAlleles, FSAReference fsaReference);
	}
}
