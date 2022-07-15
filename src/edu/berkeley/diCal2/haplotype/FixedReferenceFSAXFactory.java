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

import edu.berkeley.diCal2.haplotype.FSAXHaplotype.FSAXHapFactory;
import edu.berkeley.diCal2.haplotype.FSAXTypeConfig.FSAReference;


public class FixedReferenceFSAXFactory<H extends FSAXHaplotype> implements FSAHapFactory<H>{

	private final FSAXHapFactory<H> xFactory;
	private final FSAReference fsaRef;
	
	public FixedReferenceFSAXFactory(FSAXHapFactory<H> xFactory, FSAReference fsaRef) {
		super();
		this.xFactory = xFactory;
		this.fsaRef = fsaRef;
	}

	public H getHaplotype(int[] alleleConfig) {
		return xFactory.getHaplotype(this.fsaRef, alleleConfig);
	}

	public int getNumSegSites() {
		return this.fsaRef.getNumSegLoci();
	}
	
}
