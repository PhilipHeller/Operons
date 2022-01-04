/*
 
Copyright 2022 Philip Heller.
 
This file is part of Philip Heller's operon analysis application.

The operon analysis application is free software: you can redistribute it and/or modify it under the terms of the 
GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or 
(at your option) any later version.

The operon analysis application is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License along with the code. If not, see <https://www.gnu.org/licenses/>.
 
*/

package operons;

import java.util.*;

class Locus 
{
	String							id;
	int								start;
	int								end;
	boolean							complementStrand;
	int								contig;					// arbitrary
	String 							product;
	LinkedHashMap<Double, Double>	timeToExpressionOriginal;
	LinkedHashMap<Double, Double>	timeToExpressionNormalized;

	
	void normalize()
	{
		assert timeToExpressionOriginal != null;
		assert timeToExpressionNormalized == null;
		
		timeToExpressionNormalized = new LinkedHashMap<>();
		
		int size = timeToExpressionOriginal.size();
		double meanExpression = timeToExpressionOriginal.values().stream().mapToDouble(Double::doubleValue).sum() / size;
		for (Double time: timeToExpressionOriginal.keySet())
			timeToExpressionNormalized.put(time, timeToExpressionOriginal.get(time) - meanExpression);
	}
	
	
	boolean sameStrandAndContig(Locus that)
	{
		return this.complementStrand == that.complementStrand  &&  this.contig == that.contig;
	}
	
	
	boolean isDiel()
	{
		if (timeToExpressionOriginal == null)
			return false;
		
		double min = Collections.min(timeToExpressionOriginal.values());
		double max = Collections.max(timeToExpressionOriginal.values());
		return max - min >= 1;
	}
}
