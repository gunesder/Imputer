package edu.berkeley.path.imputer;

public class FDParameters {

	// nominal values
	private static double nom_vf = 65;			// [mph]
	private static double nom_w = 15; 			// [mph]
	private static double nom_q_max = 2000;		// [veh/mile/lane]

	private double vf;
	private double w;
	private double q_max;
	
	public FDParameters(){
		vf = FDParameters.nom_vf;
		w  = FDParameters.nom_w;
		q_max = FDParameters.nom_q_max;
	}
	
	public FDParameters(double vf,double w,double q_max){
		this.vf = vf;
		this.w  = w;
		this.q_max = q_max;
	}
	
	public void setFD(double vf,double w,double q_max){
	if(!Double.isNaN(vf))
		this.vf = vf;
	if(!Double.isNaN(w))
		this.w = w;
	if(!Double.isNaN(q_max))
		this.q_max = q_max;
	}
	
	public double getVf() {
		return vf;
	}
	
	public double getW() {
		return w;
	}
	
	public double getQ_max() {
		return q_max;
	}
	
	public double getRho_crit() {
		return q_max/vf;
	}
	
	public double getRho_jam() {
		return q_max*(1/vf+1/w);
	}
	
}
