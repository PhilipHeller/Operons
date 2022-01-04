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

import java.io.*;
import java.text.DecimalFormat;
import java.util.*;
import java.util.stream.Collectors;

import weka.core.Instance;
import weka.core.Instances;


public class StatisticsBundle extends TreeMap<Statistic, Double>
{
	public StatisticsBundle()	{ }
	
	
	public StatisticsBundle(double min, double mean, double sd, double max)
	{
		put(Statistic.MIN, min);
		put(Statistic.MEAN, mean);
		put(Statistic.SD, sd);
		put(Statistic.MAX, max);
	}
	
	
	public StatisticsBundle(Collection<Double> ds)
	{		
		double sum = ds.stream().mapToDouble(Double::valueOf).sum();
		double mean = sum / ds.size();
		double sd = ds.stream().map(d -> Math.pow(d-mean, 2)).mapToDouble(Double::valueOf).sum();
		if (sd == 0)
			sd = 0.0001;
		else
		{
			sd /= (ds.size() - 1);
			sd = Math.sqrt(sd);
		}

		put(Statistic.MIN, Collections.min(ds));
		put(Statistic.MEAN, mean);
		put(Statistic.SD, sd);
		put(Statistic.MAX, Collections.max(ds));
	}
	
	
	// For classifying an instance derived from all pairwise distances within an operon. Classifiers
	// need an Instance, but programmatic creation of an instance in poorly documented. I know how to
	// write an arff file and read it back.	
	public void toArff(File arff) throws IOException
	{
		try (FileWriter fw = new FileWriter(arff))
		{
			// Header.
			fw.write("@relation single_operon\n");
			for (Statistic s: Statistic.values())
				fw.write("@attribute " + s + " numeric\n");
			fw.write("@attribute isoperon {yes, no}\n\n@data\n");
			
			// Write a single line.
			fw.write(values().stream().map(val -> ""+val).collect(Collectors.joining(",")));
			fw.write(",yes"); 		// Hope this doesn't matter
			fw.write("\n");
		}
	}
	
	
	// Building a WEKA Instance directly is poorly documented. Play it safe, via a temporary file.
	public Instance toWEKAInstance() throws IOException
	{
		File tempf = new File(Experiment.getTempDir(), "arfftemp.arff");
		try
		{
			toArff(tempf);
			Instances i = Experiment.readArff(tempf);
			tempf.delete();
			return i.firstInstance();
		}
		catch (IOException x)
		{
			System.out.println("IO Stress: " + x.getMessage());
			return null;
		}
	}
	
	
	@Override
	public String toString()
	{
		DecimalFormat f = new DecimalFormat("#.##");
		return keySet().stream().map(s -> s + "=" + f.format(get(s))).collect(Collectors.joining("  "));
	}
	

	public String toStringForArff(String label)
	{
		return values().stream().map(val -> ""+val).collect(Collectors.joining(",")) + "," + label;
	}
	
	
	public double getMean()
	{
		return get(Statistic.MEAN);
	}
	
	
	public double getStdDev()
	{
		return get(Statistic.SD);
	}
}
