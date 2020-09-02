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

package testfrozensoil;


import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.net.URISyntaxException;
import java.util.*;
import org.jgrasstools.gears.io.timedependent.OmsTimeSeriesIteratorReader;

import bufferWriter.FreezingThawingBuffer1D;
import frozensoilsolver1D.*;
import monodimensionalProblemTimeDependent.ReadNetCDFFreezingThawing1D;
import monodimensionalProblemTimeDependent.ReadNetCDFFreezingThawingOutput1D;
//import monodimensionalProblemTimeDependent.WriteNetCDFFreezingThawing1DFloat;
import monodimensionalProblemTimeDependent.WriteNetCDFFreezingThawing1DDouble;

import org.junit.Test;

/**
 * Test the {@link Test} module.
 * To test the numerical solution against the Neumann analytical solution for 
 * the case of water freezing.
 * 
 * @author Niccolo' Tubini
 */
public class TestNeumannAnalyticalSolutionFreezing {

	@Test
	public void Test() throws Exception {


		String startDate = "2000-01-01 00:00";
		String endDate = "2000-04-10 00:00"; 
		int timeStepMinutes = 60;
		String fId = "ID";
		
		
		String pathTopBC = "resources/input/Timeseries/Neumann_268K.csv";
		String pathBottomBC = "resources/input/Timeseries/Neumann_278K.csv"; 

		String pathSaveDates = "resources/input/Timeseries/Neumann_save_all.csv"; 
		
		String pathGrid =  "resources/input/Grid_NetCDF/NeumannAnalytical_Freezing_dz_01.nc";
		String pathOutput = "resources/output/NeumannAnalytical_Freezing_dz01_3600s.nc";
	

		String outputDescription = "Solution for water freezing to compare with the Neumann analytical solution. Temperatures -5 +5, epsilon for h(T) 0.0001 [C], constant dz=0.01 [m].";
		String topBC = "Top Dirichlet";
		String bottomBC = "Bottom Dirichlet";
		String timeDelta = "3600 s";
		String sfccModel = "none";
		String stateEquationModel = "water";
		String soilThermalConductivityModel = "water";
		String interfaceThermalConductivityModel = "Harmonic mean";
		int writeFrequency = 1000;
		
		OmsTimeSeriesIteratorReader topBCReader = getTimeseriesReader(pathTopBC, fId, startDate, endDate, timeStepMinutes);
		OmsTimeSeriesIteratorReader bottomBCReader = getTimeseriesReader(pathBottomBC, fId, startDate, endDate, timeStepMinutes);

		OmsTimeSeriesIteratorReader saveDatesReader = getTimeseriesReader(pathSaveDates, fId, startDate, endDate, timeStepMinutes);


		FreezingThawingBuffer1D buffer = new FreezingThawingBuffer1D();
		WriteNetCDFFreezingThawing1DDouble writeNetCDF = new WriteNetCDFFreezingThawing1DDouble();
		ReadNetCDFFreezingThawing1D readNetCDF = new ReadNetCDFFreezingThawing1D();
		
		CallFreezingThawingSolver solver = new CallFreezingThawingSolver();
		
		
		readNetCDF.gridFilename = pathGrid;
		
		readNetCDF.read();
		
		
		solver.z = readNetCDF.z;
		solver.spaceDeltaZ = readNetCDF.spaceDelta;
		solver.volumesSoil = readNetCDF.controlVolumeDimension;
		solver.KKMAX = readNetCDF.KMAX;
		solver.temperatureIC = readNetCDF.temperatureIC;
		solver.inRheologyID = readNetCDF.rheologyID;
		solver.inParameterID = readNetCDF.parameterID;
		solver.waterDensity = 1000.0;
		solver.iceDensity = 970.0;
		solver.specificThermalCapacityWater = 4187.0;
		solver.specificThermalCapacityIce = 2108.0;
		solver.thermalConductivityWater = 0.6;
		solver.thermalConductivityIce = 2.09;
		solver.latentHeatFusion = 333700.0;
		solver.soilParticlesDensity = readNetCDF.soilParticlesDensity;
		solver.thermalConductivitySoilParticles = readNetCDF.soilParticlesThermalConductivity;
		solver.specificThermalCapacitySoilParticles = readNetCDF.soilParticlesSpecificHeatCapacity;
		solver.thetaS = readNetCDF.thetaS;
		solver.thetaR = readNetCDF.thetaR;
		solver.meltingTemperature = readNetCDF.meltingTemperature;
		solver.par1 = readNetCDF.par1; 
		solver.par2 = readNetCDF.par2; 
		solver.par3 = readNetCDF.par3;
		solver.par4 = readNetCDF.par4;
		solver.temperatureRef = 273.15;
		solver.sfccModel = sfccModel;
		solver.stateEquationModel = stateEquationModel;
		solver.soilThermalConductivityModel = soilThermalConductivityModel;
		solver.interfaceConductivityModel = interfaceThermalConductivityModel;
		solver.topBCType = topBC;
		solver.bottomBCType = bottomBC;
		solver.tTimestep = 3600;
		solver.timeDelta = 3600;
		solver.newtonTolerance = 0.003337000000000;
		solver.nestedNewton = 1;
		solver.picardIteration = 1;
		
		buffer.writeFrequency = writeFrequency;

		
		writeNetCDF.fileName = pathOutput;
		writeNetCDF.briefDescritpion = outputDescription;
		writeNetCDF.topBC = topBC;
		writeNetCDF.bottomBC = bottomBC;
		writeNetCDF.pathTopBC = pathTopBC;
		writeNetCDF.pathBottomBC = pathBottomBC;
		writeNetCDF.pathGrid = pathGrid;
		writeNetCDF.timeDelta = timeDelta;
		writeNetCDF.sfccModel = sfccModel;
		writeNetCDF.stateEquationModel = stateEquationModel;
		writeNetCDF.soilThermalConductivityModel = soilThermalConductivityModel;
		writeNetCDF.interfaceThermalConductivityModel = interfaceThermalConductivityModel;
		writeNetCDF.timeUnits= "Minutes since 01/01/1970 00:00:00 UTC";
//		writeNetCDF.mySpatialCoordinate = readNetCDF.eta; 
		writeNetCDF.mySpatialCoordinate = readNetCDF.z;
		writeNetCDF.myControlVolume = readNetCDF.controlVolumeDimension;	
		writeNetCDF.writeFrequency = writeFrequency;
		
		while( topBCReader.doProcess  ) {
		
			
			topBCReader.nextRecord();	
			HashMap<Integer, double[]> bCValueMap = topBCReader.outData;
			solver.inTopBC= bCValueMap;


			bottomBCReader.nextRecord();
			bCValueMap = bottomBCReader.outData;
			solver.inBottomBC = bCValueMap;
									
			saveDatesReader.nextRecord();
			bCValueMap = saveDatesReader.outData;
			solver.inSaveDate = bCValueMap;
			
			solver.inCurrentDate = topBCReader.tCurrent;
			
			solver.solve();

			
			buffer.inputDate = solver.inCurrentDate;
			buffer.doProcessBuffer = solver.doProcessBuffer;
			buffer.inputVariable = solver.outputToBuffer;
			
			buffer.solve();
			
			writeNetCDF.myVariables = buffer.myVariable;
			writeNetCDF.doProcess = topBCReader.doProcess;
			writeNetCDF.writeNetCDF();
			
			

		}

		topBCReader.close();
		bottomBCReader.close();
		
		/*
		 * ASSERT 
		 */
		ReadNetCDFFreezingThawingOutput1D readTestData = new ReadNetCDFFreezingThawingOutput1D();
		readTestData.filename = "resources/Output/NeumannAnalytical_Freezing_dz01_3600s_original.nc";
		readTestData.read();
		
		ReadNetCDFFreezingThawingOutput1D readSimData = new ReadNetCDFFreezingThawingOutput1D();
		readSimData.filename = pathOutput;
		readSimData.read();
	
		assertEquals(readSimData.temperature.length, readTestData.temperature.length);
		assertTrue("Temperature mismatch", Arrays.deepEquals(readSimData.temperature,readTestData.temperature));
		
	
	}

	private OmsTimeSeriesIteratorReader getTimeseriesReader( String inPath, String id, String startDate, String endDate,
			int timeStepMinutes ) throws URISyntaxException {
		OmsTimeSeriesIteratorReader reader = new OmsTimeSeriesIteratorReader();
		reader.file = inPath;
		reader.idfield = "ID";
		reader.tStart = startDate;
		reader.tTimestep = timeStepMinutes;
		reader.tEnd = endDate;
		reader.fileNovalue = "-9999";
		reader.initProcess();
		return reader;
	}
}
