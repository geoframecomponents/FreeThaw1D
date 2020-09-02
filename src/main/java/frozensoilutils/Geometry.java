package frozensoilutils;

import java.util.ArrayList;
import java.util.List;

public class Geometry {
	
	private static Geometry uniqueInstance;
	
	public static Geometry getInstance() {
		/*if (uniqueInstance == null) {
			uniqueInstance = new Variables(waterSuction, temperature);
		}*/
		return uniqueInstance;
	}
	
	public static Geometry getInstance(double[] z, double[] spaceDelta, double[] deltaZ, double surfaceElevation) {
		if (uniqueInstance == null) {
			uniqueInstance = new Geometry(z, spaceDelta, deltaZ, surfaceElevation);
		}
		return uniqueInstance;
	}
	
	
	public double[] z;
	public double[] spaceDelta;
	public double[] deltaZ;
	public double surfaceElevation;
	
	private Geometry(double[] z, double[] spaceDelta, double[] deltaZ, double surfaceElevation) {
		
		this.z = z.clone();
		this.spaceDelta = spaceDelta.clone();
		this.deltaZ = deltaZ.clone();
		this.surfaceElevation = surfaceElevation;
	}


}
