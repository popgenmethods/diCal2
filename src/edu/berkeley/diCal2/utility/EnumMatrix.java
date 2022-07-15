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

package edu.berkeley.diCal2.utility;

import java.util.EnumMap;
import java.util.Set;
import java.util.Map.Entry;

public class EnumMatrix<E extends Enum<E>, T> {
	private final EnumMap<E, EnumMap<E, T>> enumMatrix;
	private final E[] enumConstants;
	private final Class<E> keyType;
	
	public EnumMatrix(T[][] data, Class<E> keyType) {
		this.keyType = keyType;
		this.enumConstants = keyType.getEnumConstants();
		
		enumMatrix = new EnumMap<E, EnumMap<E,T>>(keyType);
		for (int rowIndex = 0; rowIndex < enumConstants.length; rowIndex++) {
			EnumMap<E,T> rowMap = new EnumMap<E,T>(keyType);
			
			for (int colIndex = 0; colIndex < enumConstants.length; colIndex++) {
				rowMap.put(enumConstants[colIndex], data[rowIndex][colIndex]);
			}
			
			enumMatrix.put(enumConstants[rowIndex], rowMap);
		}
	}
	
	public T getValue(E row, E col) {
		return enumMatrix.get(row).get(col);
	}
	
	public EnumMap<E, T> getRow(E row) {
		return enumMatrix.get(row);
	}
	
	public Class<E> getKeyType() {
		return keyType;
	}
	
	public Set<Entry<E, EnumMap<E, T>>>rowEntrySet() {
		return enumMatrix.entrySet();
	}
}
