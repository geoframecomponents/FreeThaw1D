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

package frozensoilrheology;

import closureequation.ClosureEquation;
import rheology.Rheology;

public class ClosureEquationFactory {
	
	
	public ClosureEquation create(String model, Rheology rheologyModel) {
		
		ClosureEquation myModel = null;
		
		if(model.equalsIgnoreCase("Johansen")) {
			myModel = new ThermalCondJohansen(rheologyModel);
		} else if(model.equalsIgnoreCase("Water") || model.equalsIgnoreCase("Pure Water")) {
			myModel = new ThermalCondPureWater(rheologyModel);
		} else if(model.equalsIgnoreCase("Lunardini") || model.equalsIgnoreCase("Lunardini")) {
			myModel = new ThermalCondLunardini(rheologyModel);
		} else if(model.equalsIgnoreCase("SUTRA") || model.equalsIgnoreCase("SUTRA")) {
			myModel = new SUTRAThermalCond(rheologyModel);
		} else {
			System.out.println("\n\n\tERROR: please check the soilThermalConductivityModel");
		}
		
		return myModel;
	}

}
