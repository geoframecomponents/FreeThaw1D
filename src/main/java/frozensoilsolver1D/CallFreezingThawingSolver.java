/*
 * GNU GPL v3 License
 *
 * Copyright 2019 Niccolo` Tubini
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

package frozensoilsolver1D;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.concurrent.ForkJoinPool;

import closureequation.ClosureEquation;
import frozensoilrheology.ClosureEquationFactory;
import frozensoilrheology.InternalEnergySoil;
import frozensoilrheology.StateEquationFactory;
import frozensoilutils.Bisection;
import frozensoilutils.BisectionOLD;
import frozensoilutils.Geometry;
import frozensoilutils.ProblemQuantities;
//import excessicerheology.ThermalConductivityExcessIce;
//import excessicerheology.ThermalConductivityExcessIceFactory;
import interfaceConductivity.InterfaceConductivity;
import interfaceConductivity.SimpleInterfaceConductivityFactory;
import oms3.annotations.*;
import rheology.Rheology;
//import physicalquantities.Variables;
//import soilparameters.SoilParameters;
//import rehology.RehologyParameters;
import rheology.RheologyParameters;
import rheology.SoilFreezingCharacteristicCurveFactory;
import stateequation.StateEquation;



@Description("Freezing-thawing 1D model.")
@Author(name = "Niccolo' Tubini, Stephan Gruber, and Riccardo Rigon", contact = "tubini.niccolo@gmail.com")
@Keywords("Energy equation, frozen soil, phase change, numerical solver, finite volume ")
@Documentation("")
@Bibliography("Casulli (2010)")
//@Label()
//@Name()
//@Status()
@License("General Public License Version 3 (GPLv3)")
public class CallFreezingThawingSolver {

	@Description("Number of control volumes")
	@In 
	@Unit ("-")
	public int KKMAX = 0;

	// SOIL PARAMETERS
	@Description("Water density. Default value 1000.0 [kg m-3].")
	@In 
	@Unit ("kg m-3")
	public double waterDensity = 1000.0;

	@Description("Ice density. Default value 920.0 [kg m-3].")
	@In 
	@Unit ("kg m-3")
	public double iceDensity = 920.0;

	@Description("Specific thermal capacity of water. Default value 4188.0 [J kg-1 K-1].")
	@In 
	@Unit ("J kg-1 K-1")
	public double specificThermalCapacityWater = 4188.0;

	@Description("Specific thermal capacity of ice. Default value 2117.0 [J kg-1 K-1].")
	@In 
	@Unit ("J kg-1 K-1")
	public double specificThermalCapacityIce = 2117.0;

	@Description("Thermal conductivity of water. Default value 0.6 [W m-1 K-1].")
	@In 
	@Unit ("W m-1 K-1")
	public double thermalConductivityWater = 0.6;

	@Description("Thermal conductivity of ice. Default value 2.29 [W m-1 K-1].")
	@In 
	@Unit ("W m-1 K-1")
	public double thermalConductivityIce = 2.29;

	@Description("Latent heat of fusion. Default value 333700 [J kg-1].")
	@In 
	@Unit ("J kg-1")
	public double latentHeatFusion = 333700;

	@Description("Saturated water content")
	@In 
	@Unit ("-")
	public double[] thetaS;

	@Description("Residual water content")
	@In 
	@Unit ("-")
	public double[] thetaR;

	@Description("Soil particles density")
	@In 
	@Unit ("kg m-3")
	public double[] soilParticlesDensity;

	@Description("Specific thermal capacity of soil particles")
	@In 
	@Unit ("J kg-1 K-1")
	public double[] specificThermalCapacitySoilParticles;

	@Description("Thermal conductivity of soil particles")
	@In 
	@Unit ("W m-1 K-1")
	public double[] thermalConductivitySoilParticles;

	@Description("Melting temperature")
	@In 
	@Unit ("K")
	public double[] meltingTemperature;

	@Description("First parameter of SWRC")
	@In 
	@Unit ("-")
	public double[] par1;

	@Description("Second parameter of SWRC")
	@In 
	@Unit ("-")
	public double[] par2;

	@Description("Third parameter of SWRC")
	@In 
	@Unit ("-")
	public double[] par3;

	@Description("Fourth parameter of SWRC")
	@In 
	@Unit ("-")
	public double[] par4;

	@Description("Reference temperature")
	@In 
	@Unit ("K")
	public double temperatureRef = 273.15;

	@Description("It is possibile to chose between ???? different models to compute "
			+ "the SFCC: Van Genuchten; Brooks and Corey; Kosugi unimodal")
	@In 
	public String sfccModel;

	@Description("Soil thermal conductivity model")
	@In 
	public String soilThermalConductivityModel;

	@Description("")
	@In 
	public String stateEquationModel;


	@Description("The conductivity at control volume interface can be evaluated as"
			+ " the arithmetic mean"
			+ " the maximum between"
			+ " the minimum between"
			+ " the geometric mean"
			+ " the harmonic mean")
	@In
	public String interfaceConductivityModel;

	@Description("Number of Picard iteration to update the diffusive flux matrix")
	@In
	public int picardIteration=1;
	/////////////////////////////////////////////

	@Description("Initial condition for temperature read from grid NetCDF file")
	@In
	@Unit("m")
	public double[] temperatureIC;

	@Description("z coordinate read from grid NetCDF file")
	@In
	@Unit("m")
	public double[] z;

	@Description("Space delta to compute gradients read from grid NetCDF file")
	@In 
	@Unit("m")
	public double[] spaceDeltaZ;

	@Description("Soil volume read from grid NetCDF file")
	@In 
	@Unit("m")
	public double[] volumesSoil;

	@Description("Control volume label defining the rheology")
	@In 
	@Unit("-")
	public int[] inRheologyID;

	@Description("Control volume label defining the set of the paramters")
	@In 
	@Unit("-")
	public int[] inParameterID;

	@Description("Time amount at every time-loop")
	@In
	@Unit ("s")
	public double tTimestep;

	@Description("Time step of integration")
	@In
	@Unit ("s")
	public double timeDelta;

	@Description("Tolerance for Newton iteration")
	@In
	public double newtonTolerance;

	@Description("Control parameter for nested Newton algorithm:"
			+"0 --> Newton-Raphson algorithm and the globally convergent Newton algorithm (Dall'Amico et al. 2011). The choice depends on delta."
			+"1 --> nested Newton algorithm (Casulli and Zanolli 2010)."
			+ "Default value is 1, therby nested Newton algorithm.")
	@In
	public int nestedNewton=1; 

	@Description("Control parameter of the globally convergent Newton algorithm (Dall'Amico et al. 2011). If delta=1 it is the classical "
			+ "Newton-Raphson algorithm. Default value is 1, therby Newton-Raphson algorithm")
	@In
	public double delta = 1.0; 
	// BOUNDARY CONDITIONS

	@Description("The HashMap with the time series of the boundary condition at the top of soil column")
	@In
	@Unit ("m")
	public HashMap<Integer, double[]> inTopBC;

	@Description("It is possibile to chose between 2 different kind "
			+ "of boundary condition at the top of the domain: "
			+ "- Dirichlet boundary condition --> Top Dirichlet"
			+ "- Neumann boundary condition --> Top Neumann")
	@In 
	public String topBCType;

	@Description("It is possibile to chose among 3 different kind "
			+ "of boundary condition at the bottom of the domain: "
			+ "- Dirichlet boundary condition --> Bottom Dirichlet"
			+ "- Neumann boundary condition --> Bottom Neumann")
	@In 
	public String bottomBCType;

	@Description("The HashMap with the time series of the boundary condition at the bottom of soil column")
	@In
	@Unit ("-")
	public HashMap<Integer, double[]> inBottomBC;

	@Description("The HashMap with the date to save")
	@In 
	@Unit ("")
	public  HashMap<Integer, double[]> inSaveDate;

	@Description("The current date of the simulation.")
	@In
	@Out
	public String inCurrentDate;

	@Description("ArrayList of variable to be stored in the buffer writer")
	@Out
	public ArrayList<double[]> outputToBuffer;

	@Description("ArrayList of variable to be stored in the buffer writer")
	@Out
	public boolean doProcessBuffer;


	//////////////////////////////////////////
	//////////////////////////////////////////

	@Description("Maximun number of Newton iterations")
	final int MAXITER_NEWT = 20;

	@Description("Top boundary condition according with topBCType")
	@Unit ("")
	double topBC;

	@Description("Bottom boundary condition according with bottomBCType")
	@Unit ("")
	double bottomBC;

	@Description("Newton coefficient")
	@Unit ("")
	double newtonCoeff;

	@Description("It is needed to iterate on the date")
	int step;
	int KMAX;
	double volume1;
	double volume2;
	double saveDate;

	int[] rheologyID;
	int[] parameterID;

	///////////////////////////////

	PDE1DSolver solver;
	ProblemQuantities variables;
	Geometry geometry;
	RheologyParameters rehologyParameters;
	Rheology soilFreezingCharacteristicCurve;
	SoilFreezingCharacteristicCurveFactory soilFreezingCharacteristicCurveFactory;

	@Description("This list contains the objects that describes the state equations of the problem")
	List<StateEquation> stateEquation;

	@Description("Object dealing with the internal energy model")
	StateEquation internalEnergy;
	StateEquationFactory stateEquationFactory;


	@Description("Object dealing with the thermal conductivity model")
	ClosureEquation thermalConductivity;
	ClosureEquationFactory closureEquationFactory;


	@Description("This object compute the interface hydraulic conductivity accordingly with the prescribed method.")
	InterfaceConductivity interfaceConductivity;
	SimpleInterfaceConductivityFactory interfaceConductivityFactory;

	BisectionOLD bisection;

	@Execute
	public void solve() {


		if(step==0) {

			variables = ProblemQuantities.getInstance(temperatureIC, volumesSoil, null);

			geometry = Geometry.getInstance(z, spaceDeltaZ, volumesSoil, -999.0);

			rehologyParameters = RheologyParameters.getInstance(waterDensity, iceDensity, specificThermalCapacityWater,
					specificThermalCapacityIce, thermalConductivityWater, thermalConductivityIce, latentHeatFusion, temperatureRef,
					thetaS, thetaR, soilParticlesDensity, specificThermalCapacitySoilParticles, thermalConductivitySoilParticles,
					meltingTemperature, par1, par2, par3, par4);

			outputToBuffer= new ArrayList<double[]>();

			soilFreezingCharacteristicCurveFactory = new SoilFreezingCharacteristicCurveFactory();
			soilFreezingCharacteristicCurve = soilFreezingCharacteristicCurveFactory.create(sfccModel);

			stateEquationFactory = new StateEquationFactory();
			internalEnergy = stateEquationFactory.create(stateEquationModel, soilFreezingCharacteristicCurve);

			stateEquation = new ArrayList<StateEquation>();
			stateEquation.add(internalEnergy);

			closureEquationFactory = new ClosureEquationFactory();
			thermalConductivity = closureEquationFactory.create(soilThermalConductivityModel, soilFreezingCharacteristicCurve);

			interfaceConductivityFactory = new SimpleInterfaceConductivityFactory();
			interfaceConductivity = interfaceConductivityFactory.createInterfaceConductivity(interfaceConductivityModel);

			bisection = new BisectionOLD(10e-14, inRheologyID, inParameterID, stateEquation);
			rheologyID = inRheologyID.clone();
			parameterID = inParameterID.clone();
			

			for(int k=0;k<KKMAX; k++) {
				variables.tStar[k] = bisection.findZero(273.13, 273.15 - 1e-14, k);
			}
			
			solver = new PDE1DSolver( topBCType,
					bottomBCType, interfaceConductivityModel, KKMAX, nestedNewton, newtonTolerance, delta,
					MAXITER_NEWT, picardIteration, stateEquation);
			KMAX = KKMAX;


		} // close step==0

//		System.out.println("\n"+inCurrentDate);


		doProcessBuffer = false;

		/*
		 * Get boundary conditions
		 */
		topBC = 0.0;
		topBC = (inTopBC.get(0)[0]);
		if(topBCType.equalsIgnoreCase("Top Neumann") || topBCType.equalsIgnoreCase("TopNeumann")) {
			topBC = topBC/tTimestep;
		}

		bottomBC = 0.0;
		bottomBC = inBottomBC.get(0)[0];
		if(bottomBCType.equalsIgnoreCase("Bottom Neumann") || bottomBCType.equalsIgnoreCase("BottomNeumann")) {
			bottomBC = bottomBC/tTimestep;
		}
		
		saveDate = -1.0;
		saveDate = inSaveDate.get(0)[0];

		outputToBuffer.clear();

		double sumTimeDelta = 0;

		while(sumTimeDelta < tTimestep) {
			/*
			 * volendo posso mettere qui un time-step restriction
			 * for the CN method. Remember to update lambdas
			 * for the moment semi-implicit
			 */
			if(sumTimeDelta + timeDelta>tTimestep) {
				timeDelta = tTimestep - sumTimeDelta;
			}

			sumTimeDelta = sumTimeDelta + timeDelta;


			/*
			 * Solve PDE
			 */


			variables.internalEnergy = 0.0;
			for(int i = 0; i < KMAX; i++) {
				//variables.internalEnergys[i] = stateEquation.get(rheologyID[i]).stateEquation(variables.temperatures[i], parametersID[i], i);
				variables.internalEnergys[i] = stateEquation.get(0).stateEquation(variables.temperatures[i], parameterID[i], i);
				variables.internalEnergy += variables.internalEnergys[i];
			}

			for(int picard=0; picard<picardIteration; picard++) {

				/*
				 * Compute for each control volume the thermal conductivity
				 * 
				 */

				for(int i = 0; i < KMAX; i++) {
					variables.lambdas[i] = thermalConductivity.k(variables.temperatures[i], parameterID[i], i);
					if(i==0) {
						variables.lambdasInterface[i] = variables.lambdas[i];
					}
					else {
						volume1 = variables.soilVolumes[i-1];
						volume2 = variables.soilVolumes[i];
						variables.lambdasInterface[i] = interfaceConductivity.compute(variables.lambdas[i-1],variables.lambdas[i],volume1, volume2);
					}
				}
				
				variables.lambdasInterface[KMAX] = variables.lambdas[KMAX-1];
				
				
				boolean checkData= true;
				if(checkData == false) {

					System.out.println("Thermal conductivity: ");

					for(int i=0; i<variables.lambdas.length; i++) {
						System.out.println("\t"+variables.lambdas[i]);
					}
					
					System.out.println("internal energy: ");

					for(int i=0; i<variables.internalEnergys.length; i++) {
						System.out.println("\t"+variables.internalEnergys[i]);
					}


				}
				solver.solve(topBC, bottomBC, inCurrentDate, timeDelta, KMAX, rheologyID, parameterID);
			}

			for(int i=0; i<KMAX; i++) {
				variables.thetaW[i] = soilFreezingCharacteristicCurve.f(variables.temperatures[i], parameterID[i]);
				variables.thetaI[i] = rehologyParameters.thetaS[parameterID[i]] - soilFreezingCharacteristicCurve.f(variables.temperatures[i], parameterID[i]);
			}

			/*
			 * Compute the total internal energy at time level n+1
			 * Internal energy is the sum of soil internal energy and ice
			 */
			variables.internalEnergyNew = 0.0;
			for(int k = 0; k < KMAX; k++) {		
				variables.internalEnergys[k] = stateEquation.get(0).stateEquation(variables.temperatures[k], parameterID[k], k);
				variables.internalEnergyNew += stateEquation.get(0).stateEquation(variables.temperatures[k], parameterID[k], k);
			}

			/*
			 * Compute bottom flux
			 */
			if(bottomBCType.equalsIgnoreCase("Bottom Dirichlet") || bottomBCType.equalsIgnoreCase("BottomDirichlet")){
				variables.heatFluxBottom = variables.lambdas[0]*(variables.temperatures[0]-bottomBC)/geometry.spaceDelta[0];
			} else if (bottomBCType.equalsIgnoreCase("Bottom Neumann") || bottomBCType.equalsIgnoreCase("BottomNeumann")) {
				variables.heatFluxBottom = -bottomBC;
			}

			/*
			 * Compute top flux
			 */
			if(topBCType.equalsIgnoreCase("Top Dirichlet") || topBCType.equalsIgnoreCase("TopDirichlet")){
				variables.heatFluxTop = variables.lambdas[KMAX-1]*(topBC-variables.temperatures[KMAX-1])/geometry.spaceDelta[KMAX];
			} else if (topBCType.equalsIgnoreCase("Top Neumann") || topBCType.equalsIgnoreCase("TopNeumann")) {
				variables.heatFluxTop = topBC;
			} 
			variables.errorEnergy = variables.internalEnergyNew - variables.internalEnergy - timeDelta*(variables.heatFluxTop - variables.heatFluxBottom);
//						System.out.println("\tERROR ENERGY: " + variables.errorEnergy);

		}

		/*
		 * Save temporarily outputs in a buffer
		 */
		if(saveDate == 1) {
			//			System.out.println("SaveDate = " + saveDate);

			outputToBuffer.add(variables.temperatures);
			outputToBuffer.add(variables.thetaW);
			outputToBuffer.add(variables.thetaI);
			outputToBuffer.add(variables.internalEnergys);
			outputToBuffer.add(new double[] {variables.errorEnergy});
			outputToBuffer.add(new double[] {variables.heatFluxTop});			
			outputToBuffer.add(new double[] {variables.heatFluxBottom});

			doProcessBuffer = true;
		} else {
			//			System.out.println("SaveDate = " + saveDate);
		}

		step++;
		//		System.out.println("CallEnergySolver step = " + step);

	} //// MAIN CYCLE END ////

}  /// CLOSE  ///



