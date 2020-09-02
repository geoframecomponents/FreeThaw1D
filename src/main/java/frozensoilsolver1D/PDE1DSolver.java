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
import java.util.List;
import java.util.concurrent.ForkJoinPool;

import closureequation.ClosureEquation;
import stateequation.*;
import frozensoilboundaryconditions.*;
import frozensoilrheology.*;
import frozensoilutils.*;
import interfaceConductivity.*;
import oms3.annotations.*;

@Description("This code solve the engergy equation considering the phase transition of water, latent heat of fusion, and the heat transfer is modeled using the Fourier's law."
		+ "A semi-implicit finite volume method is used to discretize the equation, and the non-linear system is solved using the nested Newton algorithm.")
@Documentation("")
@Author(name = "Niccolo' Tubini, Stephan Gruber, and Riccardo Rigon", contact = "tubini.niccolo@gmail.com")
@Keywords("Energy equation, frozen soil, phase change, numerical solver, finite volume ")
@Bibliography("")
//@Label()
//@Name()
//@Status()
@License("General Public License Version 3 (GPLv3)")
public class PDE1DSolver {


	@Description("Number of Picard iteration to update the diffusive flux matrix")
	public int picardIteration = 1;

	@Description("Time step of integration")
	@Unit ("s")
	public double timeDelta;

	@Description("Tolerance for Newton iteration")
	public double newtonTolerance = 1e-9;

	@Description("Control parameter for nested Newton algorithm:"
			+"0 --> simple Newton method"
			+"1 --> nested Newton method")
	public int nestedNewton = 1; 


	//	// BOUNDARY CONDITIONS

	@Description("It is possibile to chose between 3 different kind "
			+ "of boundary condition at the top of the domain: "
			+ "- Dirichlet boundary condition --> Top Dirichlet"
			+ "- Neumann boundary condition --> Top Neumann"
			+ "- Newton's law for heat transfer boundary condition --> Top Newton")
	public String topBCType;

	@Description("It is possibile to chose among 2 different kind "
			+ "of boundary condition at the bottom of the domain: "
			+ "- Dirichlet boundary condition --> Bottom Dirichlet"
			+ "- Neumann boundary condition --> Bottom Neumann")
	public String bottomBCType;

	//////////////////////////////////////////
	//////////////////////////////////////////

	@Description("Maximun number of Newton iterations")
	final int MAXITER_NEWT = 50;

	@Description("Top boundary condition according with topBCType")
	@Unit ("")
	double topBC;

	@Description("Bottom boundary condition according with bottomBCType")
	@Unit ("")
	double bottomBC;

	
	@Description("Vector collecting the lower diagonal entries of the coefficient matrix")
	@Unit ("?")
	double[] lowerDiagonal;

	@Description("Vector collecting the main diagonal entries of the coefficient matrix")
	@Unit ("?")
	double[] mainDiagonal;

	@Description("Vector collecting the upper diagonal entries of the coefficient matrix")
	@Unit ("?")
	double[] upperDiagonal;

	@Description("Right hand side vector of the scalar equation to solve")
	@Unit ("-")
	double[] rhss;

	@Description("Thermal conductivity at the cell interface i+1/2")
	@Unit ("W/(Km)")
	double kP;

	@Description("Thermal conductivity at the cell interface i-1/2")
	@Unit ("W/(Km)")
	double kM;

	@Description("Number of control volume for domain discetrization")
	@Unit (" ")
	int KMAX;

	//	@Description("Control variable for the integration time loop ")
	//	@Unit ("s")
	//	public double sumTimeDelta = 0.0;

	@Description("Space step")
	@Unit ("m")
	double[] spaceDelta;

	boolean checkData = false;

	@Description("Object to perform the nested Newton algortithm")
	NestedNewton nestedNewtonAlg;

	@Description("This list contains the objects that describes the state equations of the problem")
	List<StateEquation> stateEquation;

	@Description("This object compute the diagonal and right hand side entries for the uppermost cell accordingly with the prescribed top boundary condition.")
	BoundaryCondition topBoundaryCondition;

	@Description("This object compute the diagonal and right hand side entries for the lowermost cell accordingly with the prescribed bottom boundary condition.")
	BoundaryCondition bottomBoundaryCondition;

	@Description("This object compute the interface hydraulic conductivity accordingly with the prescribed method.")
	InterfaceConductivity interfaceConductivity;
	SimpleInterfaceConductivityFactory interfaceConductivityFactory;
	
	ProblemQuantities variables;
	Geometry geometry;

    //////////////////////////////



	public PDE1DSolver( String topBCType, String bottomBCType, String interfaceConductivityModel, int VECTOR_LENGTH, int nestedNewton, double newtonTolerance, double delta,
			int MAXITER_NEWT, int picardIteration, List<StateEquation>  stateEquation) {


		SimpleBoundaryConditionFactory boundCondFactory = new SimpleBoundaryConditionFactory();
		this.topBCType = topBCType;
		this.bottomBCType = bottomBCType;
		topBoundaryCondition = boundCondFactory.createBoundaryCondition(topBCType);		
		bottomBoundaryCondition = boundCondFactory.createBoundaryCondition(bottomBCType);	

		interfaceConductivityFactory = new SimpleInterfaceConductivityFactory();
		interfaceConductivity = interfaceConductivityFactory.createInterfaceConductivity(interfaceConductivityModel);

		nestedNewtonAlg = new NestedNewton(nestedNewton, newtonTolerance, MAXITER_NEWT, VECTOR_LENGTH, stateEquation, delta);


		lowerDiagonal = new double[VECTOR_LENGTH];
		mainDiagonal  = new double[VECTOR_LENGTH];
		upperDiagonal = new double[VECTOR_LENGTH];
		rhss 		  = new double[VECTOR_LENGTH];
		kP 		      = 0.0;
		kM	          = 0.0;

		variables = ProblemQuantities.getInstance();
		geometry = Geometry.getInstance();
		
	}


	/**
	 * 
	 * @param topBC
	 * @param bottomBC
	 * @param inCurrentDate
	 * @param timeDelta
	 */
	public void solve(double topBC, double bottomBC, String inCurrentDate, double timeDelta, int KMAX, int[] rheologyID, 
						int[] parameterID) {


		this.timeDelta = timeDelta;


		/* COEFFICIENT MATRIX IS BUILD BY THREE VECTORS COLLECTING ELEMENTS OF THE THREE DIAGONAL:
				   	 a lower diagonal T_(i+1)
				   	 b main diagonal  T_i
				   	 c upper diagonal T_(i-1)
				   	 RIGHT HAND SIDE 
		*/
		for(int i = 0; i < KMAX; i++) {
			if( i == 0 ) {
				
				kP = variables.lambdasInterface[i+1];
				kM = variables.lambdasInterface[i];
				lowerDiagonal[i] =  bottomBoundaryCondition.lowerDiagonal(-999.0, kP, kM, geometry.spaceDelta[i+1], geometry.spaceDelta[i], timeDelta);
				mainDiagonal[i] = bottomBoundaryCondition.mainDiagonal(-999.0, kP, kM, geometry.spaceDelta[i+1], geometry.spaceDelta[i], timeDelta);
				upperDiagonal[i] = bottomBoundaryCondition.upperDiagonal(-999.0, kP, kM, geometry.spaceDelta[i+1], geometry.spaceDelta[i], timeDelta);
				rhss[i] = variables.internalEnergys[i] + bottomBoundaryCondition.rightHandSide(bottomBC, kP, kM, geometry.spaceDelta[i+1], geometry.spaceDelta[i], timeDelta) + variables.heatSource[i];

			} else if(i == KMAX -1) {

				kP = variables.lambdasInterface[i+1];
				kM = variables.lambdasInterface[i];
				lowerDiagonal[i] = topBoundaryCondition.lowerDiagonal(-999.0, kP, kM, geometry.spaceDelta[i+1], geometry.spaceDelta[i], timeDelta); 
				mainDiagonal[i] = topBoundaryCondition.mainDiagonal(-999.0, kP, kM, geometry.spaceDelta[i+1], geometry.spaceDelta[i], timeDelta);
				upperDiagonal[i] = topBoundaryCondition.upperDiagonal(-999.0, kP, kM,  geometry.spaceDelta[i+1], geometry.spaceDelta[i], timeDelta);
				rhss[i] = variables.internalEnergys[i] + topBoundaryCondition.rightHandSide(topBC, kP, kM, geometry.spaceDelta[i+1], geometry.spaceDelta[i], timeDelta) + variables.heatSource[i];

			} else {

				kP = variables.lambdasInterface[i+1];
				kM = variables.lambdasInterface[i];
				lowerDiagonal[i] = -kM*timeDelta/geometry.spaceDelta[i];
				mainDiagonal[i] = kM*timeDelta/geometry.spaceDelta[i]  + kP*timeDelta/geometry.spaceDelta[i+1];
				upperDiagonal[i] = -kP*timeDelta/geometry.spaceDelta[i+1];
				rhss[i] = variables.internalEnergys[i] + variables.heatSource[i]; 

			}

		}
		if(checkData == true) {

			System.out.println("Lower:");
			for(int k=0;k<KMAX;k++){
				System.out.println(lowerDiagonal[k]);
			}

			System.out.println("\n\nMain:");
			for(int k=0;k<KMAX;k++){
				System.out.println(mainDiagonal[k]);
			}

			System.out.println("\n\nUpper:");
			for(int k=0;k<KMAX;k++){
				System.out.println(upperDiagonal[k]);
			}

			System.out.println("\n\nRhs:");
			for(int k=0;k<KMAX;k++){
				System.out.println(rhss[k]);
			}
		}
		/* 
		 * NESTED NEWTON ALGORITHM /
		 */
		//			s = System.currentTimeMillis();
		
		nestedNewtonAlg.set(variables.temperatures, variables.tStar, mainDiagonal, upperDiagonal, lowerDiagonal, rhss, rheologyID, parameterID, KMAX);
		variables.temperatures = nestedNewtonAlg.solver();
		
		//			stopTime = System.currentTimeMillis();
		//			elapsedTime = stopTime - startTime;
		//			System.out.println("\tNested Newton: " + elapsedTime);

		if(checkData == true) {

			System.out.println("\n\nTemperature: ");

			for(int i=0; i<KMAX; i++) {
				System.out.println("\t"+variables.temperatures[i]);
			}

		}



	} //// MAIN CYCLE END ////


}  /// CLOSE Richards1d ///



