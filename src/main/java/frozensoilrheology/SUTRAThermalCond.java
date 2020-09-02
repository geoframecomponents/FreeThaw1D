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
public class SUTRAThermalCond extends ClosureEquation{
	
	public SUTRAThermalCond(Rheology rheology) {
		super(rheology);
		// TODO Auto-generated constructor stub
	}

	private double iceRatio;
	private double kerstenNumber;
	private double lambdaSat;
	private double lambdaDry;
	
	
	public double k(double x, int id, int element) {
		
		return (1-parameters.thetaS[id])*parameters.thermalConductivitySoilParticles[id] + rheology.f(x, id)*parameters.thermalConductivityWater + (1-rheology.f(x, id))*parameters.thermalConductivityIce;
	}
	
	
}
