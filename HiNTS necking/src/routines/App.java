package routines;

import java.io.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.TimeUnit;


import classes.Sample;
import util.Configuration;
import util.Utility;

class Processor implements Callable<Double[]> {
	private int id;
	private int order_id;
	private double temp;
	private int nelec;
	private double diam;
	private int crack;
	private double screeningFactor;
	private boolean randomSeed;
	Double[] result;
	
	public Processor(int id, double temp, double diam, int nelec, int crackLength, double screeningFactor, boolean randomSeed){
		this.id = id;
		this.order_id = id;
		this.temp = temp;
		this.nelec = nelec;
		this.diam = diam;
		this.crack = crackLength;
		this.screeningFactor = screeningFactor;
		this.randomSeed = randomSeed;
		
	}
	
	public Processor(int id, int order_id, double temp, double diam, int nelec, int crackLength, double screeningFactor, boolean randomSeed){
		this.id = id;
		this.order_id = order_id;
		this.temp = temp;
		this.nelec = nelec;
		this.diam = diam;
		this.crack = crackLength;
		this.screeningFactor = screeningFactor;
		this.randomSeed = randomSeed;
	}

	@Override
	public Double[] call() {
		System.out.println("Starting: "+id);
		
		result = new Double[2];
		
		Map<String, Object> params = new HashMap<>();
		
		params.put("nelec", nelec);
		params.put("nholes", 0);
		params.put("expected_nnanops", 1200);
		params.put("sample_no", id);
		params.put("feature", "mobility");
		params.put("temperature", temp);
		params.put("closeNeighbor_thr", 0.0);
		params.put("neckedNeighbor_thr", 1.0); //0.42 allows the O_E to go up to kT
		params.put("np_diam",diam);
		params.put("crack_length", crack);
		params.put("screeningFactor", screeningFactor);
		params.put("randomSeed", randomSeed);
		
		Sample newsample = new Sample(params);
		//int copyNum = id%20;
		// the first element is the id, the second element is the actual result
		//result[0] = (double) copyNum;
		result[0] = (double) order_id;
		result[1] = newsample.simulation();
		
		System.out.println("Completed: "+id);	
		
		return result;
	}
}

public class App {

	public static void main(String[] args) throws FileNotFoundException {
		
		int numberOfSamples = 10;
		double size_std = .3;
		String stdMethod = "fixed";
		String npArraySize = "20x2x20";
		String eDensity = ".0016";
		String diameter = "6.5";
		double T = 300.0;
		boolean randomSeed = false; //whether our MersenneTwisterFast will always use the same seeds, or random seeds
		
		boolean multipleDiameters = false; //cycle through batches of NPs at different diameters
		boolean repeatCalculations = false; //repeat mobility calculations for the same samples. Use random number as seed
		boolean linearElectronFilling = false; //cycle through electron filling for samples
		//if none, then will simply go through single batch of samples
		
		if(multipleDiameters) {
			List<Integer> nelecList = Arrays.asList(76, 109, 150, 199, 259, 329, 411, 506, 614, 736, 874); //.0015e-/nm^3 .5 ll
			
			for(double screeningFactor = 1.0; screeningFactor <= 1.0; screeningFactor += 1.0) {
				int j = 0;
				
				String title = "Cubic Results" + "\\" + Configuration.latticeStructure + "_" +diameter + "nm" +npArraySize + "_" + eDensity + "_" + T + "K_screeningFactor" + screeningFactor +".txt";
				PrintWriter writer = new PrintWriter(title);
				writer.println("Diameter(nm) Mobility(cm^2/Vs) Pair_STD");
				
				for(double diam = 3; diam <= 8.0; diam = diam + 0.5) {
					System.out.println("nanoparticle diameter: " + String.valueOf(diam));
					double nelec_doub = nelecList.get(j)*.0016/.0015;
					int nelec = (int) Math.round(nelec_doub);
					ExecutorService executor = Executors.newFixedThreadPool(4);
					List<Future<Double[]>> futures = new ArrayList<>();
		
		
					int crack = 0;
					for(int i=0; i<numberOfSamples; i++){
						futures.add(executor.submit(new Processor(i,T,diam,nelec,crack,screeningFactor, randomSeed)));
					}
		
					// Stop accepting new tasks. Wait for all threads to terminate.
					executor.shutdown();
		
					System.out.println("All tasks submitted.");
		
					// wait for the processes to finish. Setting a time out limit of 1 day.
					try {
						executor.awaitTermination(1, TimeUnit.DAYS);
					} catch (InterruptedException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
		
		
					System.out.println("All tasks completed.");
		
					List<Double[]> resultRaw = new ArrayList<>();
		
					for(Future<Double[]> future: futures){
						try {
							resultRaw.add(future.get());
						} catch (InterruptedException e) {
							e.printStackTrace();
						} catch (ExecutionException e) {
							e.printStackTrace();
						}
					}
		
					// Average over regular and inverted samples
					System.out.println(resultRaw);
					double[] resultProcessed = Utility.processResult(resultRaw);
		
					double total = 0.0;
					for(int i = 0; i < resultProcessed.length; i++) {
						total += resultProcessed[i];
					}
					double mean = total/resultProcessed.length;
					
					double variance = 0.0;
					for(int i=0; i< resultProcessed.length;i++) {
						variance += Math.pow((resultProcessed[i]-mean),2);
					}
					double std = Math.sqrt(variance/(resultProcessed.length-1));
					double std_mean =  std/Math.sqrt(resultProcessed.length);
		
					String toAdd = String.valueOf(diam) + " " + mean + " " + std_mean;
					writer.println(toAdd);
					j += 1;
				}
				writer.close();
			}
		}
		
		else if(linearElectronFilling) {
			double screeningFactor = 1.0;
			String title = "Cubic Results" + "\\" + Configuration.latticeStructure + "_" +diameter + "nm" +npArraySize + "_" + eDensity + "_" + T + "K_screeningFactor" + screeningFactor +".txt";
			PrintWriter writer = new PrintWriter(title);
			writer.println("SampleNumber Mobility(cm^2/Vs)");
			
			for(double fraction = 0.7; fraction <= 0.7; fraction = fraction + 0.1) {
				int nelec = (int) Math.round(fraction*800);
				int crack = 0;
				double diam = 6.5;
				System.out.println("Number Electrons: " + String.valueOf(nelec));
				ExecutorService executor = Executors.newFixedThreadPool(4);
				List<Future<Double[]>> futures = new ArrayList<>();
	
	
	
				for(int i=0; i<numberOfSamples; i++){
					futures.add(executor.submit(new Processor(i,T,diam,nelec, crack,screeningFactor, randomSeed)));
				}
	
				// Stop accepting new tasks. Wait for all threads to terminate.
				executor.shutdown();
	
				System.out.println("All tasks submitted.");
	
				// wait for the processes to finish. Setting a time out limit of 1 day.
				try {
					executor.awaitTermination(1, TimeUnit.DAYS);
				} catch (InterruptedException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
	
	
				System.out.println("All tasks completed.");
	
				List<Double[]> resultRaw = new ArrayList<>();
	
				for(Future<Double[]> future: futures){
					try {
						resultRaw.add(future.get());
					} catch (InterruptedException e) {
						e.printStackTrace();
					} catch (ExecutionException e) {
						e.printStackTrace();
					}
				}
	
				// Average over regular and inverted samples
				System.out.println(resultRaw);
				double[] resultProcessed = Utility.processResult(resultRaw);

				for(int i = 0; i < resultProcessed.length; i++) {
					String toAdd = String.valueOf(i) + " " + resultProcessed[i];
					writer.println(toAdd);
				}
			}
			writer.close();
		}
		
		//in this scenario, we will assume that we are examining individual samples, so repeat calculations several times for each
		else if(repeatCalculations) {
			int nelec = 560; //.002 e-/nm^3 triclinic
			double diam = 6.5;
			
			int repetitions = 10;
			int iStart = 1000;
			double screeningFactor = 1.0;
			String title = "Cubic Results" + "\\" + Configuration.latticeStructure + "_" +diameter + "nm" +npArraySize + "_" + eDensity + "_" + T + "K_screeningFactor" + screeningFactor +".txt";
			
			PrintWriter writer = new PrintWriter(title);
			writer.println("SampleNumber Mobility(cm^2/Vs) Pair_STD");
			
			for(int k = 0; k < 1; k++) {
				ExecutorService executor = Executors.newFixedThreadPool(4);
				List<Future<Double[]>> futures = new ArrayList<>();
	
				//submit jobs with an ordering id separate from their actual sample id
				//this will allow us to average over groups of runs on the same sample
				for(int i=0; i<numberOfSamples; i+=2){
					for(int j=0; j<repetitions; j++) {
						int nonflippedID = i*repetitions + 2*j;
						int flippedID = i*repetitions + 2*j + 1;
						futures.add(executor.submit(new Processor(i,nonflippedID, T,diam,nelec,0,screeningFactor, randomSeed)));
						futures.add(executor.submit(new Processor(i+1,flippedID, T,diam,nelec,0, screeningFactor, randomSeed)));
					}
				}
	
				// Stop accepting new tasks. Wait for all threads to terminate.
				executor.shutdown();
	
				System.out.println("All tasks submitted.");
	
				// wait for the processes to finish. Setting a time out limit of 1 day.
				try {
					executor.awaitTermination(1, TimeUnit.DAYS);
				} catch (InterruptedException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
	
	
				System.out.println("All tasks completed.");
	
				List<Double[]> resultRaw = new ArrayList<>();
	
				for(Future<Double[]> future: futures){
					try {
						resultRaw.add(future.get());
					} catch (InterruptedException e) {
						e.printStackTrace();
					} catch (ExecutionException e) {
						e.printStackTrace();
					}
				}
				
				for(int i = 0; i < resultRaw.size(); i++) {
					System.out.println("Result: " + resultRaw.get(i));
				}
	
				// Average over regular and inverted samples
				double[] resultProcessed = Utility.processResult(resultRaw);
				System.out.println("Result processed length: " + resultProcessed.length);
	
				double[] total = new double[resultProcessed.length/repetitions];
				double[] mean = new double[total.length];
				
				for(int i = 0; i < resultProcessed.length; i++) {
					System.out.println("Result: " + resultProcessed[i]);
				}
				
				for(int i = 0; i < resultProcessed.length/repetitions; i++) {
					for(int j=0; j < repetitions; j++) {
						total[i] += resultProcessed[i*repetitions + j];
						System.out.println("Total is: " + total[i]);
					}
					mean[i] = total[i]/repetitions;
					System.out.println("Mean is: " + mean[i]);
				}
				
				double[] variance = new double[mean.length];
				double[] std_mean = new double[mean.length];
				
				for(int i=0; i< mean.length;i++) {
					for(int j = 0; j<repetitions; j++) {
						variance[i] += Math.pow((resultProcessed[i*repetitions + j]-mean[i]),2);
						System.out.println("Variance is: " + variance[i]);
					}
					std_mean[i] = Math.sqrt(variance[i]/(repetitions-1));
					std_mean[i] = std_mean[i]/Math.sqrt(repetitions);
				}
				
				for(int i = 0; i < mean.length; i++) {
					double sampleResult = mean[i];
					double sampleSTD = std_mean[i];
					String toAdd = i + " " + sampleResult + " " + sampleSTD;
					writer.println(toAdd);
				}
			}
			writer.close();
		}
		
		else {
			int nelec = 560; //.002 e-/nm^3 triclinic 20x2x20
			int nhole = 560;
			double diam = 6.5;
			double screeningFactor = 1.0;
			
			String title = "Cubic Results" + "\\" + Configuration.latticeStructure + "_" +diameter + "nm" +npArraySize + "_" + eDensity + "_" + T + "K_screeningFactor" + screeningFactor +".txt";
			
			PrintWriter writer = new PrintWriter(title);
			writer.println("SampleNumber Mobility(cm^2/Vs)");
			
			
			ExecutorService executor = Executors.newFixedThreadPool(4);
			List<Future<Double[]>> futures = new ArrayList<>();


			for(int i=0; i<numberOfSamples; i++){
				futures.add(executor.submit(new Processor(i,T,diam,nelec,0, screeningFactor, randomSeed)));
			}

			// Stop accepting new tasks. Wait for all threads to terminate.
			executor.shutdown();

			System.out.println("All tasks submitted.");

			// wait for the processes to finish. Setting a time out limit of 1 day.
			try {
				executor.awaitTermination(1, TimeUnit.DAYS);
			} catch (InterruptedException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}


			System.out.println("All tasks completed.");

			List<Double[]> resultRaw = new ArrayList<>();

			for(Future<Double[]> future: futures){
				try {
					resultRaw.add(future.get());
				} catch (InterruptedException e) {
					e.printStackTrace();
				} catch (ExecutionException e) {
					e.printStackTrace();
				}
			}

			// Average over regular and inverted samples
			double[] resultProcessed = Utility.processResult(resultRaw);

			for(int i = 0; i < resultProcessed.length; i++) {
				double sampleResult = resultProcessed[i];
				String toAdd = i + " " + sampleResult;
				writer.println(toAdd);
			}
			writer.close();
		}
		
	}

}



