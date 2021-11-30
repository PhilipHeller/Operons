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
