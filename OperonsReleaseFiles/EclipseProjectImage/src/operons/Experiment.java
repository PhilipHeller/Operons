package operons;

import java.io.*;
import java.text.DecimalFormat;
import java.util.*;
import java.util.stream.Collectors;

import org.apache.commons.math3.distribution.NormalDistribution;

import weka.classifiers.Classifier;
import weka.classifiers.Evaluation;
import weka.classifiers.bayes.BayesNet;
import weka.classifiers.bayes.NaiveBayes;
import weka.classifiers.bayes.NaiveBayesMultinomial;
import weka.classifiers.bayes.NaiveBayesMultinomialText;
import weka.classifiers.bayes.NaiveBayesMultinomialUpdateable;
import weka.classifiers.bayes.NaiveBayesUpdateable;
import weka.classifiers.functions.Logistic;
import weka.classifiers.functions.MultilayerPerceptron;
import weka.classifiers.functions.SGD;
import weka.classifiers.functions.SGDText;
import weka.classifiers.functions.SMO;
import weka.classifiers.functions.SimpleLogistic;
import weka.classifiers.functions.VotedPerceptron;
import weka.classifiers.lazy.IBk;
import weka.classifiers.lazy.KStar;
import weka.classifiers.lazy.LWL;
import weka.classifiers.meta.AdaBoostM1;
import weka.classifiers.meta.AttributeSelectedClassifier;
import weka.classifiers.meta.Bagging;
import weka.classifiers.meta.CVParameterSelection;
import weka.classifiers.meta.ClassificationViaRegression;
import weka.classifiers.meta.CostSensitiveClassifier;
import weka.classifiers.meta.FilteredClassifier;
import weka.classifiers.meta.IterativeClassifierOptimizer;
import weka.classifiers.meta.LogitBoost;
import weka.classifiers.meta.MultiClassClassifier;
import weka.classifiers.meta.MultiClassClassifierUpdateable;
import weka.classifiers.meta.MultiScheme;
import weka.classifiers.meta.RandomCommittee;
import weka.classifiers.meta.RandomSubSpace;
import weka.classifiers.meta.RandomizableFilteredClassifier;
import weka.classifiers.meta.Stacking;
import weka.classifiers.meta.Vote;
import weka.classifiers.meta.WeightedInstancesHandlerWrapper;
import weka.classifiers.misc.InputMappedClassifier;
import weka.classifiers.misc.SerializedClassifier;
import weka.classifiers.rules.DecisionTable;
import weka.classifiers.rules.JRip;
import weka.classifiers.rules.OneR;
import weka.classifiers.rules.PART;
import weka.classifiers.rules.ZeroR;
import weka.classifiers.trees.DecisionStump;
import weka.classifiers.trees.HoeffdingTree;
import weka.classifiers.trees.J48;
import weka.classifiers.trees.LMT;
import weka.classifiers.trees.REPTree;
import weka.classifiers.trees.RandomForest;
import weka.classifiers.trees.RandomTree;
import weka.core.*;
import weka.core.converters.ArffLoader;


class Experiment 
{
	private final static File					TEMP_DIR	=new File("/tmp");		// change this if you're on Windows
	
	private LinkedHashMap<String, Locus>		idToLocus;
	private List<Operon>						priorPredictedOperons;
	private List<Operon>						dielPriorPredictedOperons;	
	private List<Operon>						negativeOperons;
	private List<Operon>						dielNegativeOperons;	
	

	
	
	
							
							
							/*****************************************
							 *                                       *
							 *                TOP LEVEL              *
							 *                                       *
							 *****************************************/
	
	
	

	Experiment(File genbankFile, File timepointsCSV, File timepointsDescriptorTSV, File priorOperonsFile) 
	throws IOException
	{
		// Read loci from genbank file. 
		idToLocus = readLociFromGenbank(genbankFile);
		
		// Build time-to-expression data for measured loci.
		readExpressionMeasurements(timepointsCSV, timepointsDescriptorTSV);
		
		// Read predicted operons.
		readPriorPredictedOperons(priorOperonsFile);
		dielPriorPredictedOperons = Operon.collectMeasuredDiel(priorPredictedOperons);
		
		// Build negative training set.
		buildNegativeOperons();
		dielNegativeOperons = Operon.collectMeasuredDiel(negativeOperons);
	}
	
	
	Experiment(Experiment src)
	{
		this.idToLocus = src.idToLocus;
		this.priorPredictedOperons = new ArrayList<>(src.priorPredictedOperons);
		this.dielPriorPredictedOperons = new ArrayList<>(src.dielPriorPredictedOperons);
		this.negativeOperons = new ArrayList<>(src.negativeOperons);
		this.dielNegativeOperons = new ArrayList<>(src.dielNegativeOperons);
	}
	
	
	
			
	
	
									
									
									/*****************************************
									 *                                       *
									 *                GENBANK                *
									 *                                       *
									 *****************************************/
								
									
	
	private LinkedHashMap<String, Locus> readLociFromGenbank(File genbankFile) throws IOException
	{
		// Build a list of loci.
		int contig = 0;
		Stack<Locus> loci = new Stack<>();
		try (FileReader fr = new FileReader(genbankFile); BufferedReader br = new BufferedReader(fr))
		{
			String line;
			while ((line=br.readLine()) != null)
			{
				line = line.trim();
				if (line.startsWith("LOCUS"))
					contig++;
				if (line.startsWith("gene"))
				{
					Locus locus = new Locus();
					locus.contig = contig;
					line = line.substring(4).trim();				// complement(<4990..6219)  -OR-   4990..>6219
					line = line.replace(">", "").replace("<", "");	// complement(4990..6219)  -OR-   4990..6219
					String tagLine = br.readLine();					//                          /locus_tag="CwatDRAFT_6746"
					locus.id = tagLine.substring(tagLine.indexOf('"')+1, tagLine.lastIndexOf('"'));
					if (line.contains("complement"))
					{
						locus.complementStrand = true;
						line = line.replace("complement(", "").replace(")", "");	
					}												// 4990..6219
					for (char ch: line.toCharArray())
						assert Character.isDigit(ch) || ch == '.';
					locus.start = Integer.parseInt(line.substring(0, line.indexOf('.')));
					locus.end = Integer.parseInt(line.substring(line.lastIndexOf('.')+1));
					loci.add(locus);
				}
				else if (line.contains("/product="))
				{
					line = line.substring(line.indexOf('=') + 2);
					line = line.substring(0, line.length() - 1);
					loci.peek().product = line;
				}
			}
		}
		
		// Convert list to map.
		LinkedHashMap<String, Locus> ret = new LinkedHashMap<>();
		loci.stream().forEach(loc -> ret.put(loc.id, loc));
		return ret;
	}
	
	
	
	

	
	
						
						
						/********************************************
						 *                                          *
						 *                EXPRESSION                *
						 *                                          *
						 ********************************************/
	
	
	

	void readExpressionMeasurements(File timepointsCSV, File timepointsDescriptorTSV) throws IOException
	{
		// Descriptor file format: Each TSV line describes a column in the input spreadsheet.
		// 1st field is column number (from 0) or Excel-style column name e.g. "A" or "BX".
		// 2nd field is "geneid" or a float indicating hours since start.	
		// Any line may say "SKIP<tab>n", where n is the # of lines to skip before start of data. 
		int geneIdCol = -1;
		int nLinesToSkip = 0;
		Map<Integer, Double> colToTime = new TreeMap<>();
		try (FileReader fr = new FileReader(timepointsDescriptorTSV); BufferedReader br = new BufferedReader(fr))
		{
			String line;
			while ((line=br.readLine()) != null)
			{
				line = line.trim();
				if (line.isEmpty())
					continue;
				String[] pieces = line.split("\\t");
				if (pieces[0].equalsIgnoreCase("SKIP"))
				{
					nLinesToSkip = Integer.parseInt(pieces[1]);
					continue;
				}
				int col = -1;
				try 
				{
					col = Integer.parseInt(pieces[0]);
				}
				catch (NumberFormatException x)
				{
					col = excelColumnToInt(pieces[0]);
				}
				if (pieces[1].equalsIgnoreCase("geneid"))
					geneIdCol = col;
				else
					colToTime.put(col, Double.valueOf(pieces[1]));
			}
		}
		assert geneIdCol >= 0;
		
		// Read the spreadsheet.
		// Read the spreadsheet.
		String delim = ",";
		try (FileReader fr = new FileReader(timepointsCSV); BufferedReader br = new BufferedReader(fr))
		{
			for (int n=0; n<nLinesToSkip; n++)
				br.readLine();
			String line;
			while ((line=br.readLine()) != null)
			{
				String[] pieces = line.split(delim);
				if (line.contains("\""))
					pieces = mergeForInternalCommas(pieces);
				LinkedHashMap<Double, Double> timeToRawExpression = new LinkedHashMap<>();				
				for (Integer col: colToTime.keySet())
				{
					Double time = colToTime.get(col);
					Double xpr = Double.valueOf(pieces[col]);
					timeToRawExpression.put(time, xpr);
				}
				String id = pieces[geneIdCol];
				Locus loc = idToLocus.get(id);
				loc.timeToExpressionOriginal = timeToRawExpression;
				loc.normalize();
			}
		}	
	}
	
	
	private static int excelColumnToInt(String s)
	{
		switch (s.length())
		{
			case 1:
				return s.charAt(0) - 'A';
			case 2:
				return 26 * (s.charAt(0) - 'A' + 1) + (s.charAt(1) - 'A');
			default:
				assert false;
				return -12345;		// *sigh*
		}
	}
	
	
	// Use for split csv line with a field containing comma(s). Field must be doublequote-delimited.
	// Merge e.g. "aaaa  bbbb cccc" to "aaaa,bbbb,cccc".
	private static String[] mergeForInternalCommas(String[] originalPieces)
	{
		Stack<String> stack = new Stack<>();
		boolean merging = false;
		for (String originalPiece: originalPieces)
		{
			if (!merging)
			{
				if (originalPiece.startsWith("\"")  &&  originalPiece.endsWith("\""))
					stack.push(originalPiece.substring(1, originalPiece.length()-1));
				else if (originalPiece.startsWith("\""))
				{
					merging = true;
					stack.push(originalPiece.substring(1));
				}
				else
				{
					assert !originalPiece.contains("\"");
					stack.push(originalPiece);
				}
			}
			else
			{
				assert merging;
				if (originalPiece.endsWith("\""))
				{
					merging = false;
					originalPiece = originalPiece.substring(0, originalPiece.length()-1);
				}
				String s = stack.pop();
				s += "," + originalPiece;
				stack.push(s);
			}
		}
		return (String[])stack.toArray(new String[0]);
	}
	

	
	
	
	
	
	
	
							
							/*****************************************
							 *                                       *
							 *           PREDICTED OPERONS           *
							 *                                       *
							 *****************************************/

	
	
	//
	// Arkin format files are downloaded from the Arkin lab site at 
	//			http://www.microbesonline.org/operons/OperonList.html
	// Files are csv with one header line. Each row has 2 gene ids that are predicted to be operon partners.
	// Ids are columns 2 & 3 (from zero). Column 6 ("bOp") is "TRUE" or "FALSE". Only accept "TRUE" pairs.
	//
	private void readPriorPredictedOperons(File priorOperonsFile) throws IOException
	{
		priorPredictedOperons = new ArrayList<>();
		Map<String, Operon> idToOperon = new HashMap<>();
		
		try (FileReader fr = new FileReader(priorOperonsFile); BufferedReader br = new BufferedReader(fr))
		{
			br.readLine();				// skip header
			String line = null;
			while ((line = br.readLine()) != null)
			{
				String[] fields = line.split("\\t");
				assert fields.length >= 5;
				if (!fields[6].equals("TRUE"))
					continue;
				String id1 = fields[2];
				String id2 = fields[3];
				Locus loc1 = idToLocus.get(id1);
				Locus loc2 = idToLocus.get(id2);
				assert loc1 != null  &&  loc2 != null;
				// If neither id has been seen (almost always), create a new operon. If 1 id has been seen,
				// add the other to the operon of the seen gene. If both have been seen, might need to merge
				// 2 operons.
				int nSeen = 0;
				if (idToOperon.containsKey(id1))
					nSeen++;
				if (idToOperon.containsKey(id2))
					nSeen++;
				switch (nSeen)
				{
					case 0:
						// Neither gene has been seen. Make a new operon.
						Operon newOperon = new Operon();
						newOperon.add(loc1);
						newOperon.add(loc2);
						idToOperon.put(id1, newOperon);
						idToOperon.put(id2, newOperon);
						priorPredictedOperons.add(newOperon);
						break;
					case 1:
						// Add unseen gene to operon of seen gene.
						String seenId = idToOperon.containsKey(id1)  ?  id1  :  id2;
						String unseenId = (seenId == id1)  ?  id2  :  id1;
						Locus unseenLocus = idToLocus.get(unseenId);
						assert unseenLocus != null;
						assert !idToOperon.containsKey(unseenId);
						Operon operon = idToOperon.get(seenId);
						operon.add(unseenLocus);
						idToOperon.put(unseenId, operon);
						break;
					default:
						// Both genes have been seen. If they are in the same operon, do nothing. If they
						// are in different operons, merge the operons.
						Operon operon1 = idToOperon.get(id1);
						Operon operon2 = idToOperon.get(id2);
						if (operon1 != operon2)
						{
							Operon mergedOperon = new Operon();
							mergedOperon.addAll(operon1);
							mergedOperon.addAll(operon2);
							priorPredictedOperons.remove(operon1);
							priorPredictedOperons.remove(operon2);
							for (Locus loc: mergedOperon)
								idToOperon.put(loc.id, mergedOperon);
						}
						break;
				}
			}
		}
	}
	
							
	
	
	
								
								
								/******************************************
								 *                                        *
								 *            NEGATIVE OPERONS            *
								 *                                        *
								 ******************************************/

	
	//
	// The negative training set contains short (len = 2 or 3) predicted operons that
	// are adjacent to a strand change, plus 3 or 2 genes just past the strand switch.
	//
	private void buildNegativeOperons()
	{
		negativeOperons = new ArrayList<>();
		
		// Build.
		List<String> ids = new ArrayList<>(idToLocus.keySet());
		for (Operon op: priorPredictedOperons)
		{
			if (op.size() > 3)
				continue;
			// Check upstream.
			Locus startLoc = op.getFirstLocus();
			int startIndex = ids.indexOf(startLoc.id);
			assert startIndex >= 0;
			if (startIndex > 3)
			{
				String beforeStartId = ids.get(startIndex - 1);
				Locus beforeStartLoc = idToLocus.get(beforeStartId);
				if (startLoc.complementStrand != beforeStartLoc.complementStrand)
				{
					// Crossed a boundary => create a "noperon" of length = 5.
					Operon nop = new Operon();
					int nopStart = startIndex + op.size() - 5;
					for (int i=0; i<5; i++)
						nop.add(idToLocus.get(ids.get(nopStart + i)));
					if (nop.hasAllTxs())
						negativeOperons.add(nop);
				}
			}
			// Check downstream.
			Locus endLoc = op.getLastLocus();
			int endIndex = ids.indexOf(endLoc.id);
			assert endIndex >= 1;
			if (endIndex < ids.size() - 3)
			{
				String afterEndId = ids.get(endIndex + 1);
				Locus afterEndLoc = idToLocus.get(afterEndId);
				if (endLoc.complementStrand != afterEndLoc.complementStrand)
				{
					// Crossed a boundary => create a "noperon" of length = 5.
					startIndex = ids.indexOf(op.getFirstLocus().id);
					Operon nop = new Operon();
					for (int i=0; i<5; i++)
						nop.add(idToLocus.get(ids.get(startIndex + i)));
					if (nop.hasAllTxs())
						negativeOperons.add(nop);
				}
			}
		}
		
		// Validate.
		for (Operon op: negativeOperons)
		{
			assert op.size() == 5;
			
			boolean seenT = false;
			boolean seenF = false;
			for (Locus loc: op)
				if (loc.complementStrand)
					seenT = true;
				else
					seenF = true;
			assert seenT && seenF;
			
			for (Locus loc: op)
				assert loc.timeToExpressionOriginal != null  &&  loc.timeToExpressionNormalized != null;
		}
	}
	
	
	
	
	

	
							
							
							
							
							/*****************************************
							 *                                       *
							 *                  ARFF                 *
							 *                                       *
							 *****************************************/

	
	
	
	
	
	void writeArff(File arffFile, boolean dielOnly) throws IOException
	{
		List<Operon> positiveOps = dielOnly  ?  dielPriorPredictedOperons  :  priorPredictedOperons;
		List<Operon> negativeOps = dielOnly  ?  dielNegativeOperons  :  negativeOperons;
		
		try (FileWriter fw = new FileWriter(arffFile))
		{
			fw.write("@relation is_it_an_operon\n");
			for (Statistic s: Statistic.values())
				fw.write("@attribute " + s + " numeric\n");
			fw.write("@attribute isoperon {yes, no}\n\n@data\n");
			
			// Positive (operon) examples.
			for (Operon op: positiveOps)
			{
				if (!op.hasAllTxs())
					continue;
				StatisticsBundle stats = op.computeAblimStatistics();
				assert stats != null;
				fw.write(stats.toStringForArff("yes") + "\n");
			}
			
			// Negative (not operon) examples.
			for (Operon op: negativeOps)
			{
				StatisticsBundle stats = op.computeAblimStatistics();
				assert stats != null;
				fw.write(stats.toStringForArff("no") + "\n");
			}
		}
	}
	
	
	
	
	
	
	
	
	
								
								
								/****************************************
								 *                                      *
								 *              CLASSIFIERS             *
								 *                                      *
								 ****************************************/
							
	
	
	private final static Classifier[] CLASSIFIERS =
	{
		// bayes
		new BayesNet(), new NaiveBayes(), new NaiveBayesMultinomial(), new NaiveBayesMultinomialText(), 
		new NaiveBayesUpdateable(), new NaiveBayesMultinomialUpdateable(),
		// functions
		new Logistic(), new MultilayerPerceptron(), new SGD(), new SGDText(), new SimpleLogistic(), 
		new SMO(), new VotedPerceptron(),
		// lazy
		new IBk(), new KStar(), new LWL(),
		// meta
		new AdaBoostM1(), new AttributeSelectedClassifier(), new Bagging(), new ClassificationViaRegression(), 
		new CostSensitiveClassifier(), new CVParameterSelection(), new FilteredClassifier(), 
		new IterativeClassifierOptimizer(), new LogitBoost(), new MultiClassClassifier(), new MultiClassClassifierUpdateable(),
		new MultiScheme(), new RandomCommittee(), new RandomizableFilteredClassifier(), new RandomSubSpace(), new Stacking(), 
		new Vote(), new WeightedInstancesHandlerWrapper(),
		// misc
		new InputMappedClassifier(), new SerializedClassifier(),
		// rules
		new DecisionTable(), new JRip(), new OneR(), new PART(), new ZeroR(),
		// trees
		new DecisionStump(), new HoeffdingTree(), new J48(), new LMT(), new RandomForest(), new RandomTree(), new REPTree(),
	};
	
	
	void rankClassifiers() throws IOException					
	{
		// Build WEKA Instances.
		File arff = new File(TEMP_DIR, "arff");
		writeArff(arff, true);
		Instances instances = readArff(arff);
		arff.delete();
		
		// 5-fold cross validate all classifiers.
		Map<String, Double> classifierNameToPctCorrect = new LinkedHashMap<>();
		Set<String> failedClassifierNames = new TreeSet<>();
		for (Classifier classifier: CLASSIFIERS)
		{
			String classifierName = classifier.getClass().getName();
			classifierName = classifierName.substring(1 + classifierName.indexOf('.'));
			classifierName = classifierName.substring(1 + classifierName.indexOf('.'));
			sop(classifierNameToPctCorrect.size() + ": " + new Date() + ": " + classifierName);
			try
			{
				double pctCorrect = crossValidate(classifier, instances);
				classifierNameToPctCorrect.put(classifierName, pctCorrect);
				sop(classifierName + " " + pctCorrect);
			}
			catch (Throwable x)
			{
				failedClassifierNames.add(classifierName);
				sop(".. FAIL");
			}
		}
		sop("** Done classifying **");
		
		// Sort by accuracy.
		Map<Double, List<String>> tempMap =
			classifierNameToPctCorrect.keySet()
			.stream()
			.collect(Collectors.groupingBy(name -> classifierNameToPctCorrect.get(name)));
		Map<Double, Set<String>> pctCorrectToClassifierNames = new TreeMap<>();
		for (Double d: tempMap.keySet())
			pctCorrectToClassifierNames.put(d, new TreeSet<String>(tempMap.get(d)));
		
		// Report.
		DecimalFormat formatter = new DecimalFormat("#.##");
		for (Double accuracy: pctCorrectToClassifierNames.keySet())
			sop(formatter.format(accuracy) + ": " + pctCorrectToClassifierNames.get(accuracy).stream().collect(Collectors.joining(",")));
	}
	
	
	static double crossValidate(Classifier classifier, Instances instances) throws Exception
	{
		assert classifier != null;
		
		Evaluation eval = new Evaluation(instances);
		eval.crossValidateModel(classifier, instances, 5, new Random(0));		// 5 folds
		return eval.pctCorrect();
	}
	
	
	static Instances readArff(File f) throws IOException
	{
		FileReader fr = new FileReader(f);
		BufferedReader br = new BufferedReader(fr);
		ArffLoader.ArffReader arffr = new ArffLoader.ArffReader(br);
		Instances ret = arffr.getData();
		int nClasses = ret.numAttributes();
		ret.setClassIndex(nClasses - 1);
		return ret;
	}
	
	

	
	
	


							
							
							
							
								/*****************************************
								 *                                       *
								 *                GAUSSIAN               *
								 *                                       *
								 *****************************************/
							
							
							



	//
	// Censor each training instance in turn, retrain the LMT, and evaluate
	// the censored instance. Construct one Gaussian distribution for the 
	// positive instance scores and one for the negative instance scores.
	//
	private StatisticsBundle summarizeCensoredScoreOfNegatives() throws Exception
	{
		List<Double> censoredScores = new ArrayList<>();
		for (Operon censoredOp: dielNegativeOperons)
		{
			// Censor the operon from training, train on the rest, and evaluate the censored operon.
			Experiment censoredExper = new Experiment(this);
			censoredExper.dielPriorPredictedOperons.remove(censoredOp);		// one of these
			censoredExper.dielNegativeOperons.remove(censoredOp);			// will remove something
			Classifier lmtClassifier = buildAndTrainLTMClassifier();
			StatisticsBundle stats = censoredOp.computeAblimStatistics();
			Instance inst = stats.toWEKAInstance();
			double[] distn = lmtClassifier.distributionForInstance(inst);
			censoredScores.add(distn[0]);
		}
		return new StatisticsBundle(censoredScores);
	}	
	
	
	private NormalDistribution buildDistributionForNegativeTrainingData() throws Exception
	{
		StatisticsBundle bundle = summarizeCensoredScoreOfNegatives();
		return new NormalDistribution(bundle.getMean(), bundle.getStdDev());
	}


	//
	// Suppose an accepted item is really negative. It was accepted due to high score. What is the probability
	// that a negative item has a score this high or higher?
	//
	double probFP(NormalDistribution negativeDistribution, double scoreOfAccepted)
	{
		return 1 - negativeDistribution.cumulativeProbability(scoreOfAccepted);	// P(x >= scoreOfAccepted)
	}

	
	
	
	
	
	
	
	

							/******************************************
							 *                                        *
							 *                 MERGING                *
							 *                                        *
							 ******************************************/


	
	
	//
	// A candidate is a pair of consecutive operons, on the same strand and contig, both of length 2, each
	// having at least 1 diel gene.
	//
	List<Operon[]> collectMergeCandidates()
	{
		List<Operon[]> ret = new ArrayList<>();
		
		for (int i=0; i<priorPredictedOperons.size()-1; i++)
		{
			Operon op1 = priorPredictedOperons.get(i);
			//if (op1.size() != 2)
			//	continue;
			if (!op1.hasAllTxs())
				continue;
			if (op1.nDielGenes() == 0)
				continue;
			Operon op2 = priorPredictedOperons.get(i+1);
			//if (op2.size() != 2)
			//	continue;
			if (!op2.hasAllTxs())
				continue;
			if (op2.nDielGenes() == 0)
				continue;
			Locus lastLocusOp1 = op1.get(op1.size()-1);	
			Locus firstLocusOp2 = op2.get(0);
			if (!lastLocusOp1.sameStrandAndContig(firstLocusOp2))
				continue;
			int delta = getLoci().indexOf(firstLocusOp2) - getLoci().indexOf(lastLocusOp1);
			if (delta != 1)
				continue;
			ret.add(new Operon[] { op1, op2 } );
		}
		
		return ret;
	}
	
	
	void evaluateMergeCandidates(Classifier classifier, File otsv) throws Exception
	{
		// Evaluate candidates.
		List<Operon[]> candidates = collectMergeCandidates();
		Map<Operon[], Double> candidateToScore = new HashMap<>();
		for (Operon[] priorPair: candidates)
		{
			Operon mergedOp = new Operon(priorPair[0]);
			mergedOp.addAll(priorPair[1]);
			try
			{
				StatisticsBundle stats = mergedOp.computeAblimStatistics();
				Instance inst = stats.toWEKAInstance();
				double[] distn = classifier.distributionForInstance(inst);
				candidateToScore.put(priorPair, distn[0]);
			}
			catch (Exception x)
			{
				sop("??? " + x.getMessage());
				System.exit(1);
			}
		}
		// There's a way to do this with a single stream into a specified map type, but the
		// documentation is beyond me.
		Map<Double, List<Operon[]>> temp = 
			candidateToScore.keySet().stream().collect(Collectors.groupingBy(cand -> candidateToScore.get(cand)));		
		Map<Double, List<Operon[]>> scoreToCandidates = new TreeMap<>();
		for (Double score: temp.keySet())
			scoreToCandidates.put(score, temp.get(score));
		
		// Report.
		DecimalFormat scoreFormat = new DecimalFormat("#.##");
		DecimalFormat probFormat = new DecimalFormat("0.##E0");
		NormalDistribution negativeTrainingSetDistn = buildDistributionForNegativeTrainingData();
		try (FileWriter fw = new FileWriter(otsv))
		{
			String header = "Prior 1\tPrior 2\tScore\tP(fp)";
			fw.write(header + "\n");
			sop(header);
			for (Double score: scoreToCandidates.keySet())
			{
				for (Operon[] priors: scoreToCandidates.get(score))
				{
					double pFp = probFP(negativeTrainingSetDistn, score);
					String s = priors[0] + "\t" + priors[1] + "\t" + scoreFormat.format(score) + "\t" + probFormat.format(pFp);
					s = s.replaceAll("CwatDRAFT_", "");
					fw.write(s + "\n");
					sop(s.replaceAll("\\|", ", "));
				}
			}
		}
	}
	
	
	private Classifier buildAndTrainLTMClassifier()
	{
		Classifier classifier = new LMT();
		trainClassifier(classifier);
		return classifier;
	}
	
	
	private void trainClassifier(Classifier classifier)
	{
		try
		{
			File arff = new File(TEMP_DIR, "arff");
			writeArff(arff, true);
			Instances instances = readArff(arff);
			arff.delete();
			classifier.buildClassifier(instances);	
		}
		catch (Exception x) { }
	}
	
	
	
	
	
	
	
	
								
								
								/*****************************************
								 *                                       *
								 *              MISC & MAIN              *
								 *                                       *
								 *****************************************/


	

	List<Operon> getPredictedOperons()			{ return priorPredictedOperons; }
	List<Operon> getNegativeOperons()			{ return negativeOperons; }
	List<Locus> getLoci()						{ return new ArrayList<>(idToLocus.values()); }
	Map<String, Locus> getIdToLocusMap()		{ return idToLocus; }
	void setPredictedOperons(List<Operon> ops)	{ this.priorPredictedOperons = ops; }
	void setNegativeOperons(List<Operon> ops)	{ this.negativeOperons = ops; }
	static File getTempDir()					{ return TEMP_DIR; }
	static void sop(Object x)					{ System.out.println(x); }

	
	public static void main(String[] args) throws Exception
	{
		File gbFile = new File("data/YOUR_GENBANK_FILE");
		File timepointsCsv = new File("data/YOUR_TIMEPOINTS_FILE");
		File timepointsDescripCsv = new File("data/YOUR_TIMEPOINTS_DESCRIPTOR_FILE");
		File priorOps = new File("data/YOUR_PRIOR_OPERONS_FILE");
		Experiment exper = new Experiment(gbFile, timepointsCsv, timepointsDescripCsv, priorOps);
		
		//
		// To evaluate WEKA classifiers:
		// 		exper.rankClassifiers();
		// To classify merge candidates with an LMT model, writing a report to File otsv:
		//		Classifier classifier = new LMT();
		//      exper.evaluateMergeCandidatesWithLMT(classifier, otsv);
		// To use a different classifier model, use a different constructor call above. See static array
		// CLASSIFIERS[] for names of available WEKA classes.
		//
		
		
		
	}
}
