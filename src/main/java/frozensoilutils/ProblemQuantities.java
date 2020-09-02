package frozensoilutils;

import java.util.ArrayList;
import java.util.List;

public class ProblemQuantities {
	
	private static ProblemQuantities uniqueInstance;
	
	public static ProblemQuantities getInstance() {
		/*if (uniqueInstance == null) {
			uniqueInstance = new Variables(waterSuction, temperature);
		}*/
		return uniqueInstance;
	}
	
	public static ProblemQuantities getInstance(double[] initialCondition, double[] soilVolumes, double[] iceExcessVolumes) {
		if (uniqueInstance == null) {
			uniqueInstance = new ProblemQuantities(initialCondition, soilVolumes, iceExcessVolumes);
		}
		return uniqueInstance;
	}
	
	
	public double[] temperatures;
	public double[] tStar;
    public double[] lambdas;
	public double[] lambdasInterface;
	public double[] internalEnergys;
	public double[] soilVolumes;
//	public double[] waterExcessVolumes;
	public double[] heatSource;
//	public double[] iceExcessVolumes;
	public double[] thetaW;
	public double[] thetaI;
	
//	public double waterFromNirvana;
//	public double waterToNirvana;
//	public double energyFromNirvana;
//	public double energyToNirvana;
	
	public double errorEnergy;
//	public double errorVolume;

	public double internalEnergy;
	public double internalEnergyNew;
	public double heatFluxTop;
	public double heatFluxBottom;
	

	
	private ProblemQuantities(double[] initialCondition, double[] soilVolumes, double[] iceExcessVolumes) {
		
		temperatures = initialCondition.clone();
		tStar = new double[initialCondition.length];
		lambdas = new double[initialCondition.length];
		lambdasInterface = new double[initialCondition.length+1];
		internalEnergys = new double[initialCondition.length];
		this.soilVolumes = soilVolumes.clone();
//		waterExcessVolumes = new double[initialCondition.length];
		heatSource = new double[initialCondition.length];
//		this.iceExcessVolumes = iceExcessVolumes.clone();
		thetaW = new double[initialCondition.length];
		thetaI = new double[initialCondition.length];
		
	}


}
