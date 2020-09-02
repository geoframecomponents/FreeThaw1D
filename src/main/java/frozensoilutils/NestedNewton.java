package frozensoilutils;

import java.util.List;

import stateequation.StateEquation;

/*
 * GNU GPL v3 License
 *
 * Copyright 2017  Niccolo` Tubini
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


/**
 * This class carries out the Nested-Newton algorithm
 * (A NESTED NEWTON-TYPE ALGORITHM FOR FINITE VOLUME METHODS SOLVING RICHARDS' EQUATION IN MIXED FORM, Casulli V., Zanolli P., Journal Scientific Computing, 2010)
 *  @author Niccolo' Tubini
 */

public class NestedNewton {
	private double outerResidual;
	private double innerResidual;

	private int nestedNewton;
	private int MAXITER_NEWT;
	private int KMAX;

	private double newtonTolerance;

	private double[] x;
	private double[] dx;
	private double[] xStar;
	private double[] mainDiagonal;
	private double[] upperDiagonal;
	private double[] lowerDiagonal;
	private double[] rhss;


	private double[] fs;
	private double[] fks;
	private double[] bb;
	private double[] cc;
	private double[] dis;
	private int[] rheologyID;
	private int[] parameterID;
	private double[] x_outer;

	private double delta;

	List<StateEquation> stateEquation;
	Thomas thomasAlg = new Thomas();



	public NestedNewton(int nestedNewton, double newtonTolerance, int MAXITER_NEWT, int VECTOR_LENGTH, List<StateEquation> stateEquation, double delta){

		this.nestedNewton = nestedNewton;
		this.newtonTolerance = newtonTolerance;
		this.MAXITER_NEWT = MAXITER_NEWT;
		//		this.NUM_CONTROL_VOLUMES = VECTOR_LENGTH;
		this.stateEquation = stateEquation;
		this.delta = delta;

		x = new double[VECTOR_LENGTH];
		dx = new double[VECTOR_LENGTH];
		fs = new double[VECTOR_LENGTH];
		fks = new double[VECTOR_LENGTH];
		bb = new double[VECTOR_LENGTH]; 
		cc = new double[VECTOR_LENGTH];
		dis = new double[VECTOR_LENGTH];
		x_outer = new double[VECTOR_LENGTH];
	}




	public void set(double[] x, double[] xStar, double[] mainDiagonal, double[] upperDiagonal, double[] lowerDiagonal, double[] rhss,
			int[] rheologyID, int[] parameterID, int KMAX){

		this.x = x;
		this.xStar = xStar;
		this.mainDiagonal = mainDiagonal;
		this.upperDiagonal = upperDiagonal;
		this.lowerDiagonal = lowerDiagonal;
		this.rhss = rhss;

		this.rheologyID = rheologyID;
		this.parameterID = parameterID;
		this.KMAX = KMAX;

	}



	public double[] solver(){

		// Initial guess of T
		if(nestedNewton == 0) {

		} else {
			for(int i = 0; i < KMAX; i++) {
				x[i] = Math.min(x[i], xStar[i]-1 );
			}
		}

		//// OUTER CYCLE ////
		for(int i = 0; i < MAXITER_NEWT; i++) {
			// I have to assign 0 to outerResidual otherwise I will take into account of the previous error
			outerResidual = 0.0;
			for(int j = 0; j < KMAX; j++) {
				if(j==0) {
					fs[j] = stateEquation.get(rheologyID[j]).stateEquation(x[j],parameterID[j],j) - rhss[j]  + mainDiagonal[j]*x[j] + upperDiagonal[j]*x[j+1];
					dis[j] = stateEquation.get(rheologyID[j]).dStateEquation(x[j],parameterID[j],j);
					//					System.out.println(j+" "+fs[j]);
				} else if(j==KMAX-1) {
					fs[j] = stateEquation.get(rheologyID[j]).stateEquation(x[j],parameterID[j],j) - rhss[j] + lowerDiagonal[j]*x[j-1] + mainDiagonal[j]*x[j];
					dis[j] = stateEquation.get(rheologyID[j]).dStateEquation(x[j],parameterID[j],j);
					//					System.out.println(j+" "+fs[j]);
				} else {
					fs[j] = stateEquation.get(rheologyID[j]).stateEquation(x[j],parameterID[j],j) - rhss[j] + lowerDiagonal[j]*x[j-1] + mainDiagonal[j]*x[j] + upperDiagonal[j]*x[j+1];
					dis[j] = stateEquation.get(rheologyID[j]).dStateEquation(x[j],parameterID[j],j);
					//					System.out.println(j+" "+fs[j]);
				}

				outerResidual += fs[j]*fs[j];
			}
			outerResidual = Math.pow(outerResidual,0.5);  
			//						System.out.println("   Outer iteration " + i + " with residual " +  outerResidual);
			if(outerResidual < newtonTolerance) {
				break;
			}
			if(nestedNewton == 0){
				bb = mainDiagonal.clone();
				cc = upperDiagonal.clone();
				for(int y = 0; y < KMAX; y++) {
					bb[y] += dis[y];
				}
				thomasAlg.set(cc,bb,lowerDiagonal,fs,KMAX);
				dx = thomasAlg.solver();

				//// UPDATE SOLUTION////
				for(int s = 0; s < KMAX; s++) {
					x[s] = x[s] - dx[s]*delta; //if multiply by delta you get the globally convergent newton.
				}
			}else{
				x_outer = x.clone();
				//				for(int ii = 0; ii < KMAX; ii++) {
				//					x[ii] = Math.max(x[ii], xStar[ii] );
				//				}


				//// INNER CYCLE ////
				for(int j = 0; j < MAXITER_NEWT; j++) {
					// I have to assign 0 to innerResidual otherwise I will take into account of the previous error
					innerResidual = 0.0; 
					for(int l=0; l < KMAX; l++) {
						if(l==0) {
							fks[l] = stateEquation.get(rheologyID[l]).pIntegral(x[l],parameterID[l],l) - ( stateEquation.get(rheologyID[l]).qIntegral(x_outer[l],parameterID[l],l) + stateEquation.get(rheologyID[l]).q(x_outer[l],parameterID[l],l)*(x[l] - x_outer[l]) ) - this.rhss[l] + mainDiagonal[l]*x[l] + upperDiagonal[l]*x[l+1];
							dis[l] = ( stateEquation.get(rheologyID[l]).p(x[l],parameterID[l],l) - stateEquation.get(rheologyID[l]).q(x_outer[l],parameterID[l],l) );
							//							System.out.println(l+" "+fks[l]);
							//							System.out.println(l+" "+dis[l]);
						} else if(l==KMAX-1) {
							fks[l] = stateEquation.get(rheologyID[l]).pIntegral(x[l],parameterID[j],l) - ( stateEquation.get(rheologyID[l]).qIntegral(x_outer[l],parameterID[l],l) + stateEquation.get(rheologyID[l]).q(x_outer[l],parameterID[l],l)*(x[l] - x_outer[l]) ) - this.rhss[l] + lowerDiagonal[l]*x[l-1] + mainDiagonal[l]*x[l];
							dis[l] = ( stateEquation.get(rheologyID[l]).p(x[l],parameterID[l],l) - stateEquation.get(rheologyID[l]).q(x_outer[l],parameterID[l],l) );
							//							System.out.println(l+" "+fks[l]);
							//							System.out.println(l+" "+dis[l]);
						} else {
							fks[l] = stateEquation.get(rheologyID[l]).pIntegral(x[l],parameterID[l],l) - ( stateEquation.get(rheologyID[l]).qIntegral(x_outer[l],parameterID[l],l) + stateEquation.get(rheologyID[l]).q(x_outer[l],parameterID[l],l)*(x[l] - x_outer[l]) ) - this.rhss[l]  + lowerDiagonal[l]*x[l-1] + mainDiagonal[l]*x[l] + upperDiagonal[l]*x[l+1];
							dis[l] = ( stateEquation.get(rheologyID[l]).p(x[l],parameterID[l],l) - stateEquation.get(rheologyID[l]).q(x_outer[l],parameterID[l],l) );
							//							System.out.println(l+" "+fks[l]);
							//							System.out.println(l+" "+dis[l]);
						}

						innerResidual += fks[l]*fks[l];
					}
					innerResidual = Math.pow(innerResidual,0.5);

					//										System.out.println("     -Inner iteration " + j + " with residual " +  innerResidual);    

					if(innerResidual < newtonTolerance) {
						break;
					}

					//// THOMAS ALGORITHM////
					// Attention: the main diagonal of the coefficient matrix must not change!! The same for the upper diagonal

					bb = mainDiagonal.clone();
					cc = upperDiagonal.clone();
					for(int y = 0; y < KMAX; y++) {
						bb[y] += dis[y];
					}
					thomasAlg.set(cc,bb,lowerDiagonal,fks,KMAX);
					dx = thomasAlg.solver();

					//// UPDATE solution ////
					for(int s = 0; s < KMAX; s++) {
						x[s] = x[s] - dx[s];
						//						System.out.println(dx[s]);
					}
				}
			} //// INNER CYCLE END ////
		}
		return x;
	}

}
