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

/**
 * 
 */
package frozensoilutils;

import java.util.List;
import java.util.concurrent.RecursiveAction;

import closureequation.ClosureEquation;
import rheology.Rheology;
import rheology.RheologyParameters;
import stateequation.StateEquation;

/**
 * @author Niccolo` Tubini
 *
 */

public class Bisection extends RecursiveAction {

	private int[] rheologyID;
	private int[] parameterID;
	private List<StateEquation> stateEquation;
	private double fa;
	private double fb;
	private double fc;
	private double a;
	private double b;
	private double c;
	private double tolerance;
	private int counter;
	private double[] array;
	private int from;
	private int to;
//	public int TASK_LEN;
	public final static int TASK_LEN =25; 


	public Bisection(double[] array, int from, int to, double a, double b, double tolerance, int[] rheologyID, int[] parameterID, List<StateEquation> stateEquation) {

		this.array = array;
		this.from = from;
		this.to = to;
		this.tolerance = tolerance;
		this.rheologyID = rheologyID;
		this.parameterID = parameterID;
		this.stateEquation = stateEquation;	
//		TASK_LEN = 5;
		this.a = a;
		this.b = b;

	}

	@Override
	protected void compute() {
		int len = to - from;
		if (len <= TASK_LEN) {
			findZero(array,from, to, a, b);
		} else {
			// split work in half, execute sub-tasks asynchronously
			int mid = (from + to) >>> 1;
			//			new Bisection2(array, from, mid, a, b, tolerance, rheologyID, parameterID, stateEquation).fork();
			//			new Bisection2(array, from, mid, a, b, tolerance, rheologyID, parameterID, stateEquation).fork();
			invokeAll(new Bisection(array, from, mid, a, b, tolerance, rheologyID, parameterID, stateEquation),
					new Bisection(array, mid, to, a, b, tolerance, rheologyID, parameterID, stateEquation));
		}
	}

	private void findZero(double[] array, int from, int to, double a, double b) {

		double tmp = -9999.0;
		double aa;
		double bb;
		for (int element = from; element < to; element++) {
			counter = 1;
			aa = a;
			bb = b;
//			System.out.println("\t aa:"+a+" bb:"+ b );
//			System.out.println(parameterID[element]);
			fa = stateEquation.get(rheologyID[element]).ddStateEquation(aa, parameterID[element], element);
			fb = stateEquation.get(rheologyID[element]).ddStateEquation(bb, parameterID[element], element);
			c = (aa+bb)/2;
			fc = stateEquation.get(rheologyID[element]).ddStateEquation(c, parameterID[element], element);

			if(fc == 0.0) {
				tmp = stateEquation.get(rheologyID[element]).computeXStar(parameterID[element], element);
				if(tmp == -9999.0) {
					System.out.println("\t Error computing T_star");
				}
				array[element] = tmp;
//				//				System.out.println("\t\t" +tmp);
//				System.out.println("\t\t" +(stateEquation.get(rheologyID[element]).computeXStar(parameterID[element], element)==tmp));
//				System.out.println("\t\t" + (array[element]==tmp));
//				System.out.println("\t\tarray[element] " + array[element]);

			} else {

				while(Math.abs(aa-bb) > tolerance) {

					if(fa*fb > 0.0) {
						System.out.println("\tBISECTION: error finding zero.");

					}

					c = (aa+bb)/2;
					fc = stateEquation.get(rheologyID[element]).ddStateEquation(c, parameterID[element], element);

					if(fc == 0.0) {
						array[element] = c;
					} else {
						if(fa*fc < 0.0) {
							bb = c;
							fb = fc;
						} else {
							aa = c;
							fa = fc;
						}
					}

					if(counter>250) {
						System.out.println("\tBISECTION: reached 250 iteration |a-b| = "+Math.abs(aa-bb));
						array[element] = c;
					}
					counter ++;
				}
				
				array[element] = c;
				System.out.println("element: " +element+" parameterID: " +parameterID[element]+ "counter: "+counter+", intervall: "+Math.abs(aa-bb)+ "Tstar: "+array[element]);
			}
		}
	}
}
