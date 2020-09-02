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
package frozensoilrheology;

import closureequation.ClosureEquation;
import rheology.Rheology;
import rheology.RheologyParameters;

/**
 * @author Niccolo` Tubini
 *
 */
public class ThermalCondLunardini extends ClosureEquation{
	
	public ThermalCondLunardini(Rheology rheology) {
		super(rheology);
		// TODO Auto-generated constructor stub
	}

	
	
	public double k(double x, int id, int element) {
		
		if(x>=273.15) {
			return 2.417196;
		} else if(x<=super.parameters.meltingTemperature[id])  {
			return 3.462696;
		} else {
			return 2.939946;
		}
		
	}
	
	
}
