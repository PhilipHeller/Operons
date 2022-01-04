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


public class Ablim 
{
	public static double computeAblimDistance(LinkedHashMap<Double, Double> txs1, LinkedHashMap<Double, Double> txs2)
	{
		assert txs1.keySet().equals(txs2.keySet());
		
		// Add areas.
		double area = 0;
		List<Double> timesList = new ArrayList<>(txs1.keySet());
		for (int nStart=0; nStart<txs1.size()-1; nStart++)
		{
			Double tStart = timesList.get(nStart);
			double expr0Start = txs1.get(tStart);
			double expr1Start = txs2.get(tStart);
			Double tEnd = timesList.get(nStart+1);
			double expr0End = txs1.get(tEnd);
			double expr1End = txs2.get(tEnd);
			area += areaBetweenSegments(tStart, expr0Start, expr1Start, tEnd, expr0End, expr1End);
		}
		return area;
	}
	
	
	// y0a and y0b have common x. Ditto y1a and y1b. Returns true if segment a crosses segment b.
	private static boolean segmentsIntersect(double yStartA, double yStartB, double yEndA, double yEndB)
	{
		// If segments meet without crossing at either endpoint, they don't strictly intersect.
		if (yStartA == yStartB  ||  yEndA == yEndB)
			return false;
		
		return (yStartA < yStartB) != (yEndA < yEndB);
	}
	

	// Returns coordinates of intersection of 2 lines that are defined by member points. Lines are
	// x1--x2 and x3--x4
	// Source: http://en.wikipedia.org/wiki/Line-line_intersection
	private static double[] intersect2lines(double x1, double y1, double x2, double y2,
											double x3, double y3, double x4, double y4)
	{
		double[] ret = new double[2];		// { x, y }
		
		// x
		double numer = (x1*y2 - y1*x2) * (x3-x4)  -  (x1-x2) * (x3*y4 - y3*x4);
		double denom = (x1-x2) * (y3-y4)  -  (y1-y2) * (x3-x4);
		ret[0] = numer / denom; 
		
		// y
		numer = (x1*y2 - y1*x2) * (y3-y4)  -  (y1-y2) * (x3*y4 - y3*x4);
		ret[1] = numer / denom;
		
		return ret;
	}
		
	
	// Segments A and B have same starting/ending xs.
	private static double areaBetweenSegments(double xStart, double yStartA, double yStartB, 
											  double xEnd, double yEndA, double yEndB)
	{		
		// If segments intersect, divide area into 2 triangles and recurse.
		if (segmentsIntersect(yStartA, yStartB, yEndA, yEndB))
		{
			double[] intersection = intersect2lines(xStart, yStartA, xEnd, yEndA, 
												    xStart, yStartB, xEnd, yEndB);
			double leftArea = 
				areaBetweenSegments(xStart, yStartA, yStartB, intersection[0], intersection[1], intersection[1]);
			double rightArea = 
				areaBetweenSegments(intersection[0], intersection[1], intersection[1], xEnd, yEndA, yEndB);
			return leftArea + rightArea;
		}

		// Simple case: trapezoid.
		double deltaX = xEnd - xStart;
		assert deltaX >= -.1f  :  deltaX + " at " + new Date();
		double meanHeight = ( (yStartB-yStartA) + (yEndB-yEndA)) / 2;
		meanHeight = Math.abs(meanHeight);
		return deltaX * meanHeight;
	}


}
