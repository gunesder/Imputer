package edu.berkeley.path.imputer;

import java.util.*;

public class FDCalibrator {
	
	// fields
	double vf;			// [mph]
	double w;			// [mph]
	double q_max;		// [veh/hour/lane]
	double rho_crit;	// [veh/mile/lane]
	double w_min = 8;	// [mph]
	double w_max = 19;	// [mph]
		
	// getters and setters
	
	// constructors
	
	// methods
	public Detector calibrateParameters(Detector d){
		// organize into an array of DataPoint
		ArrayList<DataPoint> datavec = new ArrayList<DataPoint>();
		int i;
		for(i=0;i<d.getDensityData().size();i++)
			datavec.add(new DataPoint(d.getDensityData().get(i), d.getFlowData().get(i), d.getSpeedData().get(i)));
		
		// maximum flow and its corresponding density
		DataPoint maxflw = new DataPoint(0,Double.NEGATIVE_INFINITY,0);
		for(i=0;i<d.getFlowData().size();i++)
			if(datavec.get(i).flw>maxflw.flw)
				maxflw.setval(datavec.get(i));

		q_max = maxflw.flw;
		
		// split data into congested and freeflow regimes ...............
		ArrayList<DataPoint> congestion = new ArrayList<DataPoint>();		// congestion states
		ArrayList<DataPoint> freeflow = new ArrayList<DataPoint>();			// freeflow states
		for(i=0;i<d.getDensityData().size();i++)
			if(datavec.get(i).dty>=maxflw.dty)
				congestion.add(datavec.get(i));
			else
				freeflow.add(datavec.get(i));
		
		// vf is the average freeflow speed
		vf = percentile("spd",freeflow,0.5f);
		
		// compute critical density
		rho_crit = q_max/vf;

		// BINNING
		ArrayList<DataPoint> supercritical = new ArrayList<DataPoint>(); 	// data points above rho_crit
		for(i=0;i<d.getDensityData().size();i++)
			if(datavec.get(i).dty>=rho_crit)
				supercritical.add(datavec.get(i));

		// sort supercritical w.r.t. density
		Collections.sort(supercritical);

		int numsupercritical = supercritical.size();
		int Bin_width = 10;
		int step=Bin_width;
		ArrayList<DataPoint> BinData = new ArrayList<DataPoint>();
		for(i=0;i<numsupercritical;i+=Bin_width){

			if(i+Bin_width>=numsupercritical)
		        step = numsupercritical-i;

		    if(step!=0){
		    	List<DataPoint> Bin = (List<DataPoint>) supercritical.subList(i,i+step);
		    	if(!Bin.isEmpty()){
			        double a = 2.5f*percentile("flw",Bin,0.75f) - 1.5f*percentile("flw",Bin,0.25f); 			        
			        double b = percentile("flw",Bin,1f);
			        BinData.add(new DataPoint(percentile("dty",Bin,0.5f),Math.min(a,b),Float.NaN));
		    	}
		    }
		}

		// Do constrained LS
		ArrayList<Double> ai = new ArrayList<Double>();
		ArrayList<Double> bi = new ArrayList<Double>();
		for(i=0;i<BinData.size();i++){
			bi.add(q_max - BinData.get(i).flw);
			ai.add(BinData.get(i).dty - rho_crit);
		}

		if(BinData.size()>0){
			float sumaibi = 0;
			float sumaiai = 0;
			for(i=0;i<BinData.size();i++){
				sumaibi += ai.get(i)*bi.get(i);
				sumaiai += ai.get(i)*ai.get(i);
			}
			w = sumaibi/sumaiai;
			w = Math.max(w,w_min);
			w = Math.min(w,w_max);
		}
		else{
		    w  = Float.NaN;
		}

		// store parameters in sensor
		d.setFdParams(new FDParameters((float) vf, (float) w, (float) q_max));
		return d;		
		
	}
	
	private static Double percentile(String qty,List<DataPoint> x,float p){
		// compute the p'th percentile qty (p in [0,1])
		ArrayList<Double> values = new ArrayList<Double>();
		int numdata = x.size();
		if(qty.equals("spd"))
			for(int i=0;i<numdata;i++)
				values.add(x.get(i).spd);
		if(qty.equals("flw"))
			for(int i=0;i<numdata;i++)
				values.add(x.get(i).flw);
		if(qty.equals("dty"))
			for(int i=0;i<numdata;i++)
				values.add(x.get(i).dty);

		Collections.sort(values);

		if(p==0)
			return values.get(0);
		if(p==1)
			return values.get(numdata-1);

		int z = (int) Math.floor(numdata*p);
		if(numdata*p==z)
			return (values.get(z-1)+values.get(z))/2f;
		else
			return values.get(z);
	}
	
	// nested classes
	private static class DataPoint implements Comparable<Object> {
		double dty;
		double flw;
		double spd;
		public DataPoint(double d,double f,double s){
			dty=d;
			flw=f;
			spd=s;
		}
		public void setval(DataPoint x){
			this.dty=x.dty;
			this.flw=x.flw;
			this.spd=x.spd;
		}
		public int compareTo(Object x) {
			double thatdty = ((DataPoint) x).dty;
			if(this.dty==thatdty)
				return 0;
			else if(this.dty>thatdty)
				return 1;
			else
				return -1;
		}
	}

}
