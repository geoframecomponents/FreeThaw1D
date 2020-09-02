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
public class InternalEnergySoil extends StateEquation {

	private Rheology soilModel;
	protected ProblemQuantities variables;
	private double f;
	private double df;
	private double ddf;

	public InternalEnergySoil(Rheology rheology) {
		this.soilModel = rheology;
		this.variables = ProblemQuantities.getInstance();
	}



	@Override
	public double stateEquation(double x, int id, int element) {
		
		f = soilModel.f(x, id);
		return ( super.parameters.specificThermalCapacitySoilParticles[id]*super.parameters.soilParticlesDensity[id]*(1.0-super.parameters.thetaS[id]) + 
				super.parameters.waterDensity*super.parameters.specificThermalCapacityWater*f + 
				super.parameters.iceDensity*super.parameters.specificThermalCapacityIce*(super.parameters.thetaS[id]-f) )*variables.soilVolumes[element]*(x-super.parameters.temperatureRef)
				+ super.parameters.waterDensity*super.parameters.latentHeatFusion*f*variables.soilVolumes[element];

	}


	@Override
	public double dStateEquation(double x, int id, int element) {

		f = soilModel.f(x, id);
		df = soilModel.df(x, id);
		return ( super.parameters.specificThermalCapacitySoilParticles[id]*super.parameters.soilParticlesDensity[id]*(1.0-super.parameters.thetaS[id]) + 
				super.parameters.waterDensity*super.parameters.specificThermalCapacityWater*f + 
				super.parameters.iceDensity*super.parameters.specificThermalCapacityIce*(super.parameters.thetaS[id]-f) )*variables.soilVolumes[element]
						+ (super.parameters.waterDensity*super.parameters.specificThermalCapacityWater - super.parameters.iceDensity*super.parameters.specificThermalCapacityIce)*df*variables.soilVolumes[element]*(x-super.parameters.temperatureRef)
						+ super.parameters.waterDensity*super.parameters.latentHeatFusion*df*variables.soilVolumes[element];
	}


	@Override
	public double ddStateEquation(double x, int id, int element) {
//		System.out.println("ddStateEquation: "+id);
		f = soilModel.f(x, id);
		df = soilModel.df(x, id);
		ddf = soilModel.ddf(x, id);
		return 2*(super.parameters.waterDensity*super.parameters.specificThermalCapacityWater - super.parameters.iceDensity*super.parameters.specificThermalCapacityIce)*df*variables.soilVolumes[element] +
				(super.parameters.waterDensity*super.parameters.specificThermalCapacityWater - super.parameters.iceDensity*super.parameters.specificThermalCapacityIce)*ddf*variables.soilVolumes[element]*(x-super.parameters.temperatureRef) 
				+ super.parameters.waterDensity*super.parameters.latentHeatFusion*ddf*variables.soilVolumes[element];
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
			return stateEquation(variables.tStar[element], id, element) + dStateEquation(variables.tStar[element], id, element)*(x-variables.tStar[element]);
		}

	}



	@Override
	public double computeXStar(int id, int element) {
		// TODO Auto-generated method stub
		return -9999.0;
	}



}
