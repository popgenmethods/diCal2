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

public class Pair<A extends Comparable<A>, B extends Comparable<B>> implements Comparable<Pair<A,B>>	{
	private A first;
	private B second;
	
	public Pair(A first, B second)	{	
		this.first = first;	
		this.second = second;	
	}
	
	public A first()	{	
		return first;	
	}
	
	public B second()	{	
		return second;	
	}

	public int compareTo(Pair<A, B> other) {
		int res = first.compareTo(other.first);
		if (res == 0)	res = second.compareTo(other.second);
		return res;
	}
	
	public boolean equals(Object o)	{
		if (o != null && this.getClass() == o.getClass())	{
			Pair<?, ?> other = (Pair<?, ?>) o;
			return ((first == other.first || (first != null && other.first != null && first.equals(other.first)))	&&
					(second == other.second || (second != null && other.second != null && second.equals(other.second)))
					);
		}
		return false;
	}
	
	public int hashCode() {
        int hashFirst = first != null ? first.hashCode() : 0;
        int hashSecond = second != null ? second.hashCode() : 0;

        return (hashFirst * 0x1f1f1f1f) ^ hashSecond;
    }
}
