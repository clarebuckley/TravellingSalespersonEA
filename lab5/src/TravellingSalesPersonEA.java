import java.util.ArrayList;
import java.util.Collections;
import java.util.Random;
import java.util.Scanner;
import java.io.File;
import java.io.FileNotFoundException;

public class TravellingSalesPersonEA {


	private ArrayList<ArrayList<Integer>> population;
	private final int populationSize;
	private final double mutationProbability;
	private final int numberOfGenerations;
	private final int tournamentSize;

	private double[][] csvDistances = createEuclideanGraph();

	public TravellingSalesPersonEA(int populationSizeParam, double mutationProbabilityParam, int generationsParam, int tournamentSizeParam) throws FileNotFoundException {
		populationSize = populationSizeParam;
		mutationProbability = mutationProbabilityParam;
		numberOfGenerations = generationsParam;
		tournamentSize = tournamentSizeParam;

		System.out.println("Initialising population of " + populationSize);
		System.out.println("--> mutation probability: " + mutationProbability);
		System.out.println("--> tournament size: " + tournamentSize);
		initialisePopulation();
		
		System.out.println("\nGoing through " + numberOfGenerations + " generations...");
		evolve();
		
		System.out.println("\nFinding best cost from final population...");
		System.out.println("----------------------------");
		double bestCost = evaluateFinalPopulation();
		System.out.println("Best cost: " + bestCost);
		System.out.println("----------------------------");

	}

	public static void main(String[] args) throws FileNotFoundException {
		new TravellingSalesPersonEA(200, 0.70, 6000, 100);
	}

	public void evolve() {
		int generations = 0;
		while(generations < numberOfGenerations) {
			oneGeneration();
			generations++;
		}
	}

	public double evaluateFinalPopulation() {
		//Find best route from end population
		double bestCost = 10000;
		ArrayList<Integer> bestRoute = new ArrayList<Integer>();
		for(int i = 0; i < population.size(); i++) {
			ArrayList<Integer> thisRoute = population.get(i);
			double thisCost = getCostOfRoute(thisRoute);
			if(thisCost < bestCost) {
				bestCost = thisCost;
				bestRoute = thisRoute;
			}
		}
		System.out.println("Best route: " + bestRoute); 
		return bestCost;
	}


	public void oneGeneration() {
		//Select parents
		ArrayList<Integer> parent1 = tournamentParentSelection();
		ArrayList<Integer> parent2 = tournamentParentSelection();
		//Recombine parents
		ArrayList<Integer> child = generateCrossover(parent1, parent2);
		//Mutate resulting offspring and add to possible solutions
		if(Math.random() < mutationProbability) {
			child =	generateTwoSwap(child);
		}

		ArrayList<ArrayList<Integer>> newPopulation = population;
		newPopulation.add(child);
		//Replace weakest member of population
		population = replaceWeakestIndividual(newPopulation);

	}


	/**
	 * Initialise population for evolutionary algorithm
	 */
	private void initialisePopulation() {
		population = new ArrayList<ArrayList<Integer>>();
		for(int i = 0; i < populationSize; i++) {
			population.add(generateRandomRoute());
		}
	}

	/**
	 * Select one parent from the route with lowest cost from
	 * a subgroup of the population
	 * @return
	 */
	private ArrayList<Integer> tournamentParentSelection(){
		Random random = new Random();
		int sampleSize = tournamentSize;
		ArrayList<ArrayList<Integer>> candidates = new ArrayList<ArrayList<Integer>>();

		//Get candidate parents
		for(int i = 0; i < sampleSize; i++) {
			int randomIndex = random.nextInt(populationSize);
			candidates.add(population.get(randomIndex));
		}

		//Get best candidate from selection
		ArrayList<Integer> bestCandidate = new ArrayList<Integer>();
		double bestCost = 10000;
		for(ArrayList<Integer>candidate : candidates) {
			double thisCost = getCostOfRoute(candidate);
			if(thisCost < bestCost) {
				bestCost = thisCost;
				bestCandidate = candidate;
			}
		}

		return bestCandidate;
	}


	private ArrayList<Integer> generateCrossover(ArrayList<Integer> parent1, ArrayList<Integer> parent2){
		Random random = new Random();
		ArrayList<Integer> child = new ArrayList<Integer>(parent1.size());
		int crossoverSize = (int)parent1.size()/2;
		//Set up empty child to be filled in
		for(int i = 0; i < parent1.size(); i++) {
			child.add(-1);
		}

		//Get unique indexes to be used in crossover
		ArrayList<Integer> indexes = new ArrayList<Integer>();
		while(indexes.size() < crossoverSize) {
			int randomIndex = random.nextInt(parent1.size());
			if(!indexes.contains(randomIndex)) {
				indexes.add(randomIndex);
			}
		}

		//Find values of parent 1 to be added to child
		ArrayList<Integer> tmpParent2 = new ArrayList<Integer>(parent2);
		for(int randomIndex : indexes) {
			child.set(randomIndex, parent1.get(randomIndex));
			int tmpParent2Index = tmpParent2.indexOf(parent1.get(randomIndex));
			if(tmpParent2Index >= 0) {
				tmpParent2.set(tmpParent2Index, -1);
			}
		}

		//Get values of parent 2 > 0 to be added to child
		ArrayList<Integer> addToChildren = new ArrayList<Integer>();
		for(int i = 0; i < parent2.size(); i++) {
			if(tmpParent2.get(i) >= 0) {
				addToChildren.add(parent2.get(i));
			}
		}

		//Add parent 2 values to child
		int j = 0;
		for(int i = 0; i < child.size(); i++) {
			if(child.get(i) == -1 && addToChildren.size() > 0) {
				child.set(i, addToChildren.get(j));
				j++;
			}
		}

		return child;
	}


	/**
	 * Swap two cities in a route to generate a neighbour of routeInput
	 * 
	 * @param routeInput
	 *            - route to swap cities in
	 * @return new route with two swapped cities
	 */
	private ArrayList<Integer> generateTwoSwap(ArrayList<Integer> route){
		Random random = new Random();
		int swap1 = random.nextInt(route.size());
		int swap2 = random.nextInt(route.size());
		while(swap1 == swap2){
			swap2 = random.nextInt(route.size());
		}
		Collections.swap(route, swap1, swap2);
		return route;
	}


	private ArrayList<ArrayList<Integer>> replaceWeakestIndividual(ArrayList<ArrayList<Integer>> candidates){
		double highestCost = 0;
		ArrayList<Integer> weakestCandidate = new ArrayList<Integer>();
		for(ArrayList<Integer> candidate : candidates) {
			double thisCost = getCostOfRoute(candidate);
			if(getCostOfRoute(candidate) > highestCost) {
				highestCost = thisCost;
				weakestCandidate = candidate;
			}
		}
		int weakIndex = population.indexOf(weakestCandidate);
		candidates.remove(weakIndex);
		return candidates;
	}


	/**
	 * Calculate cost of route for a route generated by csv
	 * 
	 * @param route
	 *            city route to find cost for
	 * @return
	 */
	private double getCostOfRoute(ArrayList<Integer> route) {
		double totalCost = 0;
		for (int i = 0; i < route.size() - 1; i++) {
			Integer start = route.get(i);
			Integer next = route.get(i + 1);

			double val = csvDistances[start][next];
			totalCost = totalCost + val;
		}

		// Loop route back to start point
		int firstCity = route.get(0);
		int lastCity = route.get(route.size() - 1);
		totalCost = totalCost + csvDistances[lastCity][firstCity];
		return totalCost;
	}

	/**
	 * Generates new legal random route based on data from a given CSV file
	 * 
	 * @return String containing generated route
	 */
	private ArrayList<Integer> generateRandomRoute() {
		ArrayList<Integer> cities = new ArrayList<Integer>();

		for(int i = 0; i < csvDistances.length; i++){
			cities.add(i);
		}

		Collections.shuffle(cities);
		return cities;
	}



	/**
	 * Parse the city values from ulysses.csv and extract distance values for
	 * new 2d array
	 * 
	 * @param file
	 * @return 2d array of x/y values for each city location
	 * @throws FileNotFoundException
	 */
	private double[][] parseCsv(String file) throws FileNotFoundException {
		Scanner scanner = new Scanner(new File(file));
		scanner.useDelimiter("\n");
		double[][] csvGraph = new double[16][2];

		while (scanner.hasNext()) {
			char[] line = scanner.next().trim().toCharArray();
			if (Character.isDigit(line[0])) {
				String[] lineElements = new String(line).split(",");
				int index = Integer.parseInt(lineElements[0]) - 1;
				csvGraph[index][0] = Double.parseDouble(lineElements[1]);
				csvGraph[index][1] = Double.parseDouble(lineElements[2]);
			}
		}
		scanner.close();
		return csvGraph;
	}

	/**
	 * Find the Euclidean distance between cities A and B
	 * 
	 * @param ax
	 *            x value of 'from' city
	 * @param ay
	 *            y value of 'from' city
	 * @param bx
	 *            x value of 'to' city
	 * @param by
	 *            y value of 'to' city
	 * @return - euclidean distance from city A to city B
	 */
	private double euclideanDistance(double ax, double ay, double bx, double by) {
		return Math.sqrt(Math.pow(bx - ax, 2) + (Math.pow(by - ay, 2)));
	}

	/**
	 * Create 16x16 graph to show the euclidean distance between a route given
	 * in a csv file
	 * 
	 * @return
	 * @throws FileNotFoundException
	 */
	private double[][] createEuclideanGraph() throws FileNotFoundException {
		double[][] euclideanGraph = new double[16][16];
		double[][] csvGraph = parseCsv(
				"C:\\Users\\Clare\\Documents\\Aston\\Y3\\Computational Intelligence\\Lab 1\\src\\ulysses16.csv");

		for (int x = 0; x < 16; x++) {
			for (int y = 0; y < 16; y++) {
				// x or y relate to a and b in euclideanDistance, 0 and 1 relate
				// to x and y in euclideanDistance
				euclideanGraph[x][y] = euclideanDistance(csvGraph[x][0], csvGraph[x][1], csvGraph[y][0],
						csvGraph[y][1]);
			}
		}
		return euclideanGraph;
	}


	// Method to find factorial of given number
	static int factorial(int n) {
		int res = 1, i;
		for (i = 2; i <= n; i++) {
			res *= i;
		}
		return res;
	}
}
