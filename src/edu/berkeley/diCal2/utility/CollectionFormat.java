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

import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.NavigableMap;
import java.util.TreeMap;
import java.util.Map.Entry;

public class CollectionFormat {
	public static String formatDoubleCollection(Collection<Double> collection, DecimalFormat format, String sep, String start, String end) {
		StringBuffer sb = new StringBuffer(start);

		for (double item : collection) {
			sb.append(format.format(item) + sep);
		}
		
		if (collection.size() > 0) {
			sb.replace(sb.length() - sep.length(), sb.length(), end);
		} else {
			sb.append(end);
		}
		
		return sb.toString();
	}
	
	public interface Formatter<T> {
		public String formatItem(T item);
	}
	
	public static <T> String formatCollection(Collection<T> collection, String sep, String start, String end) {
		return formatCollection(collection, sep, start, end, null);
	}
	
	public static <T> String formatCollection(Collection<T> collection, String sep, String start, String end, Formatter<T> formatter) {
		StringBuffer sb = new StringBuffer(start);

		for (T item : collection) {
			sb.append((formatter == null ? item : formatter.formatItem(item)) + sep);
		}
		
		if (collection.size() > 0) {
			sb.replace(sb.length() - sep.length(), sb.length(), end);
		} else {
			sb.append(end);
		}
		
		return sb.toString();
	}
	
	public static <T> String formatArray(T[] array, String sep, String start, String end) {
		return formatArray(array, sep, start, end, null);
	}
	
	public static <T> String formatArray(T[] array, String sep, String start, String end, Formatter<T> formatter) {
		StringBuffer sb = new StringBuffer(start);

		for (T item : array) {
			sb.append((formatter == null ? item : formatter.formatItem(item)) + sep);
		}
		
		if (array.length > 0) {
			sb.replace(sb.length() - sep.length(), sb.length(), end);
		} else {
			sb.append(end);
		}
		
		return sb.toString();
	}
	
	
	
	public static String formatArray(double[] array, String sep, String start, String end) {
		StringBuffer sb = new StringBuffer(start);

		for (double item : array) {
			sb.append(item + sep);
		}
		
		if (array.length > 0) {
			sb.replace(sb.length() - sep.length(), sb.length(), end);
		} else {
			sb.append(end);
		}
		
		return sb.toString();
	}

	public static String formatArray(double[] array, String sep, String start, String end, DecimalFormat format) {
		StringBuffer sb = new StringBuffer(start);
		
		for (double item : array) {
			sb.append(format.format(item) + sep);
		}
		
		if (array.length > 0) {
			sb.replace(sb.length() - sep.length(), sb.length(), end);
		} else {
			sb.append(end);
		}
		
		return sb.toString();
	}
	
	public static String formatArray(int[] array, String sep, String start, String end) {
		StringBuffer sb = new StringBuffer(start);

		for (int item : array) {
			sb.append(item + sep);
		}
		
		if (array.length > 0) {
			sb.replace(sb.length() - sep.length(), sb.length(), end);
		} else {
			sb.append(end);
		}
		
		return sb.toString();
	}
	
	public static String formatMapArray(int[] array, String sep, String mapsTo, String start, String end) {
		NavigableMap<Integer, List<Integer>> cMap = new TreeMap<Integer, List<Integer>>();
		
		for (int aIndex = 0; aIndex < array.length; aIndex++) {
			List<Integer> siteList = cMap.get(array[aIndex]);
			
			if (siteList == null) {
				siteList = new ArrayList<Integer>();
				cMap.put(array[aIndex], siteList);
			}
			
			siteList.add(aIndex);
		}
		
		StringBuffer sb = new StringBuffer(start);
		
		for (Entry<Integer, List<Integer>> item : cMap.entrySet()) {
			sb.append(formatArray(item.getValue().toArray(), " ", "", "") + mapsTo + item.getKey() + sep);
		}
		
		if (cMap.size() > 0) {
			sb.replace(sb.length() - sep.length(), sb.length(), end);
		} else {
			sb.append(end);
		}
		
		return sb.toString();
	}

	public static<E,F> String formatMap(Map<E, F> map, String sep, String mapsTo, String start, String end) {
		StringBuffer sb = new StringBuffer(start);
		
		for (Entry<E, F> item : map.entrySet()) {
			sb.append(item.getKey() + mapsTo + item.getValue() + sep);
		}

		if (map.size() > 0) {
			sb.replace(sb.length() - sep.length(), sb.length(), end);
		} else {
			sb.append(end);
		}
		
		return sb.toString();
	}
	
	public static<E,F> String formatMap(Map<E, F> map, String sep, String mapsTo, String start, String end, Formatter<F> formatVal) {
		StringBuffer sb = new StringBuffer(start);
		
		for (Entry<E, F> item : map.entrySet()) {
			sb.append(item.getKey() + mapsTo + formatVal.formatItem(item.getValue()) + sep);
		}

		if (map.size() > 0) {
			sb.replace(sb.length() - sep.length(), sb.length(), end);
		} else {
			sb.append(end);
		}
		
		return sb.toString();
	}
	
	public static<E,F> String formatRMap(Map<E, F> map, String sep, String mapsTo, String start, String end) {
		Map<F, List<E>> rMap = new HashMap<F, List<E>>();
		
		for (Entry<E, F> mapEntry : map.entrySet()) {
			List<E> siteList = rMap.get(mapEntry.getValue());
			
			if (siteList == null) {
				siteList = new ArrayList<E>();
				rMap.put(mapEntry.getValue(), siteList);
			}
			
			siteList.add(mapEntry.getKey());
		}
		
		StringBuffer sb = new StringBuffer(start);
		
		for (Entry<F, List<E>> item : rMap.entrySet()) {
			sb.append(formatArray(item.getValue().toArray(), " ", "", "") + mapsTo + item.getKey() + sep);
		}
		
		if (rMap.size() > 0) {
			sb.replace(sb.length() - sep.length(), sb.length(), end);
		} else {
			sb.append(end);
		}
		
		return sb.toString();
	}
}
