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

import frozensoilutils.ProblemQuantities;
import rheology.Rheology;
import stateequation.*;

/**
 * @author Niccolo` Tubini
 *
 */
public class InternalEnergyLunardini extends StateEquation {

	private Rheology soilModel;
	protected ProblemQuantities variables;

	public InternalEnergyLunardini(Rheology rheology) {
		this.soilModel = rheology;
		this.variables = ProblemQuantities.getInstance();
	}



	@Override
	public double stateEquation(double x, int id, int element) {
		
		return 690030*(x-273.15)*variables.soilVolumes[element] + 
				super.parameters.iceDensity*super.parameters.latentHeatFusion*soilModel.f(x,id)*variables.soilVolumes[element];
	}


	@Override
	public double dStateEquation(double x, int id, int element) {

		return 690030*variables.soilVolumes[element] + 
				super.parameters.iceDensity*super.parameters.latentHeatFusion*soilModel.df(x, id)*variables.soilVolumes[element];
	}


	@Override
	public double ddStateEquation(double x, int id, int element) {

		return 0.0;
	}


	@Override
	public double p(double x, int id, int element) {

		if(x<=variables.tStar[element]) {
			return dStateEquation(x, id, element);  
		} else {
			return dStateEquation(variables.tStar[element], id, element);
		}

	}


	@Override
	public double pIntegral(double x, int id, int element) {

		if(x<=variables.tStar[element]) {
			return stateEquation(x, id, element);  
		} else {
			return stateEquation(variables.tStar[element], id, element) + p(variables.tStar[element], id, element)*(x-variables.tStar[element]);
		}
		
	}



	@Override
	public double computeXStar(int id, int element) {
		// TODO Auto-generated method stub
		return super.parameters.meltingTemperature[id];
	}



}
