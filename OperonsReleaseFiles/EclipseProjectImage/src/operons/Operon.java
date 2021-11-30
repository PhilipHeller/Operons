package operons;

import java.util.*;
import java.util.stream.Collectors;


class Operon extends ArrayList<Locus>
{
	Operon()	{ }
	
	
	Operon(Collection<Locus> src)
	{
		super(src);
	}
	
	
	@Override
	public String toString()
	{
		return this.stream().map(loc->loc.id).collect(Collectors.joining("|"));
	}
	
	
	Locus getFirstLocus()
	{
		return get(0);
	}
	

	Locus getLastLocus()
	{
		return get(size()-1);
	}
	
	
	boolean hasAllTxs()
	{
		for (Locus loc: this)
			if (loc.timeToExpressionOriginal == null)
				return false;
		return true;
	}
	
	
	StatisticsBundle computeAblimStatistics()
	{
		List<Double> ds = new ArrayList<>();
		for (int i=0; i<size()-1; i++)
		{
			Locus loc1 = get(i);
			Locus loc2 = get(i+1);
			double d = Ablim.computeAblimDistance(loc1.timeToExpressionNormalized, loc2.timeToExpressionNormalized);
			ds.add(d);
		}
		return new StatisticsBundle(ds);
	}
	
	
	int nDielGenes()
	{
		int n = 0;
		for (Locus loc: this)
			if (loc.isDiel())
				n++;
		return n;
	}
	
	
	Operon[] split(Random rand)
	{
		if (size() < 4)
			return null;
		
		int maxLenA = size() - 2;
		int lenA = randomInRange(rand, 2, maxLenA);
		Operon opA = new Operon();
		for (int i=0; i<lenA; i++)
			opA.add(get(i));
		Operon opB = new Operon();
		for (int i=lenA; i<size(); i++)
			opB.add(get(i));
		return new Operon[] { opA, opB };
	}
	
	
	List<String> getIds()
	{
		return this.stream().map(loc -> loc.id).collect(Collectors.toList());
	}
		
	
	static int randomInRange(Random r, int min, int max)		// E.g. 5-10
	{
		int nPossibleVals = max - min + 1;						// 6
		int delta = Math.abs(r.nextInt()) % nPossibleVals;		// r.nextInt(int) seems to have problems
		return min + delta;										// 5 + [0..5] = [5..10]
	}
	
	
	static List<Operon> collectMeasuredDiel(List<Operon> src)
	{
		return src.stream()
			.filter(op -> op.hasAllTxs())
			.filter(op -> op.nDielGenes() >= 0)
			.collect(Collectors.toList());
	}
	
	
	static void sop(Object x)	{ System.out.println(x); }
	
	
	public static void main(String[] args) throws Exception
	{
		sop("START");
		
		Operon op = new Operon();
		for (int i=0; i<6; i++)
			op.add(new Locus());
		for (int seed=0; seed<20; seed++)
		{
			Random r = new Random(seed);
			Operon[] ops = op.split(r);
			sop(ops[0].size() + ":" + ops[1].size());
		}
		
		sop("DONE");
	}
}
