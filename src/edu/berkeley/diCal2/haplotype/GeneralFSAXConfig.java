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


public abstract class GeneralFSAXConfig<H extends FSAXHaplotype> extends GeneralGeneticTypeConfig<H> implements FSAXTypeConfig<H> {
	public final int numLoci;
	public final int numAlleles;
	public final FSAReference fsaReference;
	
	public GeneralFSAXConfig(int numLoci, int numAlleles, FSAReference fsaReference) {
		this.numLoci = numLoci;
		this.numAlleles = numAlleles;
		this.fsaReference = fsaReference;
	}
	
	public int numLoci() { return numLoci; }
	public int numAlleles() { return numAlleles; }
	
	public FSAReference getReference() {
		return fsaReference;
	}
	
	public final void adjustType(H type, int adjustment) {
		if (type.getReference() != fsaReference) {
			throw new UnsupportedOperationException("adjusting a type with different reference unsupported!");
		}
		
		super.adjustType(type, adjustment);
	}
}
