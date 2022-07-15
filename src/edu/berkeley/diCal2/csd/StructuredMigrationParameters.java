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

import java.util.ArrayList;

import edu.berkeley.diCal2.csd.auxiliary.MyEigenDecomposition;
import edu.berkeley.diCal2.demography.Demography;
import edu.berkeley.diCal2.utility.RealPartition.Interval;

public class StructuredMigrationParameters {
	
	final TrunkProcess trunk;
	final private Demography demo;
	final private ArrayList<MyEigenDecomposition> eigenStuffList;

	public StructuredMigrationParameters (TrunkProcess trunk) {
		
		this.trunk = trunk;
		
		// get the demo stuff out of demography
		this.demo = trunk.getRefinedDemography();
		
		// just some testing
		assert (trunk.getSampleSizes().length == demo.numAncientDemes(0));

		// get the extended matrices
		this.eigenStuffList = new ArrayList<MyEigenDecomposition>();
		for (int i=0; i<this.demo.migrationMatrixList.size(); i++) {
			if (this.demo.isPulse(i)) {
				// add null, cause it's a pulse
				this.eigenStuffList.add (null);
			}
			else {
				// get the extended matrix Z
				double[][] tmpExtMatrix = ExtendedMigrationMatrix.computeExtendedMigrationMatrix (this.demo.migrationMatrixList.get(i), this.trunk.getAbsorbRates(i));
				
				// get the eigendecomposition
				MyEigenDecomposition tmpDecomp = new MyEigenDecomposition(tmpExtMatrix, UberDemographyCore.ONELOCUS_EPSILON);
				// check that there is no eigenvalue with zero real and non-zero imaginary part
				for (int j=0; j<tmpDecomp.lambda.length; j++) {
					assert (!((Math.abs(tmpDecomp.lambda[j].re) < UberDemographyCore.ONELOCUS_EPSILON) && (Math.abs(tmpDecomp.lambda[j].im) > UberDemographyCore.ONELOCUS_EPSILON)));
				}
				assert (tmpDecomp.lambda.length/2 == this.demo.migrationMatrixList.get(i).length);
				// should be a fine thing now
	
				// and add it to list
				this.eigenStuffList.add (tmpDecomp);
			}
		}
		
		// all done now, right?
	}
	
	public Interval[] getIntervalList() {
		return this.demo.epochList;
	}
	
	public int numIntervals () {
		return this.demo.epochList.length;
	}
	
	public ArrayList<MyEigenDecomposition> getEigenStuffList () {
		return eigenStuffList;
	}

	public int numAncientDemes (int interval) {
		assert (interval <= this.demo.numIntervals());
		return this.demo.numAncientDemes(interval);
	}

	public Demography getDemography(){
		return demo;
	}
}

