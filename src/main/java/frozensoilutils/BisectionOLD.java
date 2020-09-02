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

import closureequation.ClosureEquation;
import rheology.Rheology;
import rheology.RheologyParameters;
import stateequation.StateEquation;

/**
 * @author Niccolo` Tubini
 *
 */

public class BisectionOLD {

	private int[] rheologyID;
	private int[] parameterID;
	private List<StateEquation> stateEquation;
	private double fa;
	private double fb;
	private double fc;
	private double a;
	private double b;
	private double c;
	private double tmp;
	private double tolerance;
	private int counter;

	public BisectionOLD(double tolerance, int[] rheologyID, int[] parameterID, List<StateEquation> stateEquation) {

		this.tolerance = tolerance;
		this.rheologyID = rheologyID;
		this.parameterID = parameterID;
		this.stateEquation = stateEquation;		

	}



	public double findZero(double a, double b, int element) {

		counter = 1;

		fa = stateEquation.get(rheologyID[element]).ddStateEquation(a, parameterID[element], element);
		fb = stateEquation.get(rheologyID[element]).ddStateEquation(b, parameterID[element], element);
		c = (a+b)/2;
		fc = stateEquation.get(rheologyID[element]).ddStateEquation(c, parameterID[element], element);

		if(fc == 0.0) {
			
			tmp = stateEquation.get(rheologyID[element]).computeXStar(parameterID[element], element);
			
			if(tmp == -9999.0) {
				System.out.println("\t Error computing T_star");
			}
			
			c = stateEquation.get(rheologyID[element]).computeXStar(parameterID[element], element);

		} else {
			while(Math.abs(a-b) > tolerance) {


				if(fa*fb > 0.0) {
					System.out.println("\tBISECTION: error finding zero.");
				}

				c = (a+b)/2;
				fc = stateEquation.get(rheologyID[element]).ddStateEquation(c, parameterID[element], element);

				if(fc == 0.0) {
					return c;
				} else {
					if(fa*fc < 0.0) {
						b = c;
						fb = fc;
					} else {
						a = c;
						fa = fc;
					}
				}

				if(counter>150) {
					System.out.println("\tBISECTION: reached 250 iteration |a-b| = "+Math.abs(a-b));
					return c;
				}
				counter ++;
			}
		}
//		System.out.println("element: " +element+ "counter: "+counter+", intervall: "+Math.abs(a-b));

		return c;

	}
}
