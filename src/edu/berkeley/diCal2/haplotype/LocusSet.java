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
import java.io.OutputStream;
import java.io.PrintStream;
import java.io.Reader;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

@SuppressWarnings("serial")
public class LocusSet extends ArrayList<Integer> {
	public LocusSet(Collection<Integer> locusSet) {
		super(locusSet);
	}
	
	public static LocusSet readLocusSet(Reader reader) throws IOException {
		BufferedReader bufferedReader = new BufferedReader(reader);
		String firstLine = bufferedReader.readLine();
		
		String[] locusSetString = firstLine.split(" ");

		List<Integer> locusSet = new ArrayList<Integer>(locusSetString.length);
		for (int locusIdx = 0; locusIdx < locusSetString.length; locusIdx++) {
			locusSet.add(Integer.parseInt(locusSetString[locusIdx]));
		}

		return new LocusSet(locusSet);
	}
	
	public static void writeLocusSet(LocusSet outLocusSet, OutputStream oStream) {
		PrintStream printStream = new PrintStream(oStream);
		
		StringBuilder locusString = new StringBuilder();
		for (int locusIdx = 0; locusIdx < outLocusSet.size(); locusIdx++) {
			locusString.append(outLocusSet.get(locusIdx));
			if (locusIdx < outLocusSet.size() - 1) locusString.append(",");

			printStream.println(locusString);
		}
	}
	
	public static LocusSet limitSiteList(LocusSet siteList, int maxSites) {
		List<Integer> newSegregatingSites = new ArrayList<Integer>(maxSites);
		
		if (siteList.size() <= maxSites) {
			newSegregatingSites.addAll(siteList);
		} else {
			int skipSites = (siteList.size() - maxSites) / 2;
			for (int siteIdx = skipSites; siteIdx < skipSites + maxSites; siteIdx++) {
				newSegregatingSites.add(siteList.get(siteIdx));
			}
		}
		
		return new LocusSet(newSegregatingSites);
	}
}
