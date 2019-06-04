package classes;

import java.util.ArrayList;
import java.util.Arrays;

import util.Configuration;
import util.Constants;

public class Electron extends Charge {
	public Electron(){
		charge = -1.0*Constants.sqrt2;
		mass = Configuration.emass;
		ratesOnCharge = 0.0;
		hoppings = new ArrayList<HoppingEvent>();
		visitedNPs = new ArrayList<double[]>();
		//events = new RegularHopping[5];
	}
	
	
	
	public double lookForEvents(Sample sample) {
		//System.out.println(this.hostNP);
		double degeneracy;
		double netDistance;
		double newRate;
		HoppingEvent newEvent;
		Nanoparticle sourceNP = this.hostNP;
		int sourceOrbital = this.hostOrbital;
		
		this.hoppings.clear();
		
		// zero total rates on charge
		ratesOnCharge = 0 ;

		// add cluster hoppings next
		for(Nanoparticle neighborNP : sourceNP.closeNeighbors){
			for(int band=0; band<sample.nbnd; band++){
				// band is targetOrbital
				// allow hopping towards drain only
				netDistance = (neighborNP.z-sourceNP.z-sample.cellz*Math.round((neighborNP.z-sourceNP.z)/sample.cellz));
				if(netDistance>0 && neighborNP.occupationCB[band]<neighborNP.occupationMAX[band] ){
					
					newEvent = new HoppingEvent(this, neighborNP, sourceNP, band, sourceOrbital, true);
					degeneracy = (double) (neighborNP.occupationMAX[band] - neighborNP.occupationCB[band]) ; 
					newRate = -charge*degeneracy*(sample.metalPrefactor)*sample.voltage/sample.cellz*netDistance;
					newEvent.setRate(newRate);
					this.hoppings.add(newEvent);
					ratesOnCharge += newRate ;
					//System.out.println("CLUSTER RATE IS: " + newRate);
				}
			}
		}
		
		//add necking transport next
		if(sourceNP.neckedNeighbors.length != 0) {
			for(Nanoparticle neighborNP : sourceNP.neckedNeighbors) {
				for(int band = 0; band < sample.nbnd; band++) {
					//band is target orbital
					double energy_diff = 0.0;
					
					// Kinetic energy difference
			    	energy_diff += neighborNP.cbenergy[band] - sourceNP.cbenergy[sourceOrbital] ;
			    	
			    	// On-Site charging energy
			    	energy_diff += calculateCharging(neighborNP, sourceNP, sample);
			    	
			    	// Electron-hole exciton interaction
			    	energy_diff += calculateExciton(neighborNP, sourceNP, sample);
	
					// allow hopping only towards nanoparticles with a lower energy level
					netDistance = (neighborNP.z-sourceNP.z-sample.cellz*Math.round((neighborNP.z-sourceNP.z)/sample.cellz));
					//System.out.println(charge*sample.voltage*netDistance/sample.cellz);
					energy_diff += charge*sample.voltage*netDistance/sample.cellz;
					//energy_diff < 0
					
					double neckRadius = sourceNP.neckRadiusMap.get(neighborNP);
					double diameter = (sourceNP.diameter + neighborNP.diameter)/2.0;
					double overlap_energy = sample.ediff_neck_thr*2*neckRadius/diameter; //Scale proportionally to neck thickness
					//System.out.println("Overlap energy is: " + overlap_energy);
					//System.out.println("kT is: " + sample.temperature);
					//double overlap_energy = sample.temperature
					
					if(energy_diff <= overlap_energy && neighborNP.occupationCB[band]<neighborNP.occupationMAX[band] ){
						
						newEvent = new HoppingEvent(this, neighborNP, sourceNP, band, sourceOrbital, false);
						degeneracy = (double) (neighborNP.occupationMAX[band] - neighborNP.occupationCB[band]) ; 
						
						double thisNP_electron_density = sourceNP.occupationTotalElectron/(Math.PI/6.0*Math.pow(sourceNP.diameter, 3.0));
						double neighborNP_electron_density = (neighborNP.occupationTotalElectron + 1)/(Math.PI/6.0*Math.pow(neighborNP.diameter, 3.0));
						//double electron_density = sample.nelec/(sample.cellx*sample.celly*sample.cellz);
						double electron_density = (thisNP_electron_density + neighborNP_electron_density)/2.0;
						double matrixElement = 9*electron_density*Math.pow(neckRadius, 3.0)/(Configuration.emass*Math.pow(diameter, 2.0));  //9*hbar^2*n*rho^3/(m^* * d^2)
						
						newRate = 2*Math.PI*Math.pow(matrixElement, 2.0)*degeneracy;
						newEvent.setRate(newRate);
						this.hoppings.add(newEvent);
						ratesOnCharge += newRate ;
					}
					else if(energy_diff > overlap_energy && neighborNP.occupationCB[band]<neighborNP.occupationMAX[band]) {
						newEvent = new HoppingEvent(this, neighborNP, sourceNP, band, sourceOrbital, false);
						degeneracy = (double) (neighborNP.occupationMAX[band] - neighborNP.occupationCB[band]) ; 
						double thisNP_electron_density = sourceNP.occupationTotalElectron/(Math.PI/6.0*Math.pow(sourceNP.diameter, 3.0));
						double neighborNP_electron_density = (neighborNP.occupationTotalElectron + 1)/(Math.PI/6.0*Math.pow(neighborNP.diameter, 3.0));
						//double electron_density = sample.nelec/(sample.cellx*sample.celly*sample.cellz);
						double electron_density = (thisNP_electron_density + neighborNP_electron_density)/2.0;
						double matrixElement = 9*electron_density*Math.pow(neckRadius, 3.0)/(Configuration.emass*Math.pow(diameter, 2.0));
						
						double thermal_activation = Math.exp(-energy_diff/sample.temperature);
						//System.out.println("Thermal activation is: " + thermal_activation);
						newRate = 2*Math.PI*Math.pow(matrixElement, 2.0)*degeneracy*thermal_activation;
						newEvent.setRate(newRate);
						this.hoppings.add(newEvent);
						ratesOnCharge += newRate ;
					}
				}
			}
		}
			
		// add regular hoppings last
		for(Nanoparticle neighborNP : sourceNP.nearestNeighbors){
			if(!Arrays.asList(sourceNP.neckedNeighbors).contains(neighborNP)) {
				for(int band=0; band<sample.nbnd; band++){
					if(neighborNP.occupationCB[band]<neighborNP.occupationMAX[band]){
						//HoppingEvent(Electron e, Nanoparticle target, Nanoparticle source, int to, int so, boolean intra)
						degeneracy = (double) (neighborNP.occupationMAX[band] - neighborNP.occupationCB[band]) ;
						
						double npdistance = neighborNP.edgeDistanceMap.get(sourceNP);
						double overlap = Math.sqrt(-sample.emass*(sourceNP.cbenergy[sourceOrbital] + neighborNP.cbenergy[band]) / 2.0) ;
						//System.out.println("Hopping prefactor is: " + degeneracy*sample.jumpfreq * Math.exp(-2.0*npdistance*overlap));
						//degeneracy = 1;
						newEvent = new HoppingEvent(this, neighborNP, sourceNP, band, sourceOrbital, false);
						newRate = degeneracy*calculateRate(neighborNP, sourceNP, band, sourceOrbital, sample);
						newEvent.setRate(newRate);
						//System.out.println("Hopping rate is: " + newRate);
						this.hoppings.add(newEvent);
						ratesOnCharge += newRate;
					}
				}
			}
		}
		
		return ratesOnCharge;
	}

	
	private double calculateCharging(Nanoparticle targetNP, Nanoparticle sourceNP, Sample sample){
    	// On-Site charging energy
		/*OLD METHOD
		 * 
		 * double charging = 0.0;
	    
    	if(sample.lcapacitance0)
    		charging += targetNP.selfenergy0 - sourceNP.selfenergy0 ;
    	
    	charging += targetNP.occupationTotalElectron*targetNP.selfenergy - (sourceNP.occupationTotalElectron-1)*sourceNP.selfenergy;
    	
	    */
		
		//new method
	    double chargingBefore = (targetNP.calculateOnSiteElectronCharging(targetNP.occupationTotalElectron)
                + sourceNP.calculateOnSiteElectronCharging((sourceNP.occupationTotalElectron)));
	    double chargingAfter = (targetNP.calculateOnSiteElectronCharging(targetNP.occupationTotalElectron+1) 
	    		+ sourceNP.calculateOnSiteElectronCharging((sourceNP.occupationTotalElectron-1)));


	    return (chargingAfter-chargingBefore)/sample.screeningFactor;
	    
		//return charging/Configuration.screeningFactor;
    }
	
	private double calculateExciton(Nanoparticle targetNP, Nanoparticle sourceNP, Sample sample){
	    /*old method for calculating excitonic interactions
	     * 
	     * double exciton = 0.0;
    	
	    double excitonSP = Math.min(sourceNP.occupationTotalHoles, sourceNP.occupationTotalElectron);
	    double excitonTP = Math.min(targetNP.occupationTotalHoles, targetNP.occupationTotalElectron+1);
	    
	    exciton = (excitonTP>0) ? -targetNP.selfenergy0 : 0 ;
	    exciton += (excitonSP>0) ? sourceNP.selfenergy0 : 0 ; 
	    
	    excitonSP = Math.min(sourceNP.occupationTotalHoles, sourceNP.occupationTotalElectron-1);
	    excitonTP = Math.min(targetNP.occupationTotalHoles, targetNP.occupationTotalElectron);

	    exciton += -excitonTP*targetNP.selfenergy + excitonSP*sourceNP.selfenergy ;
	    
    	return exciton;*/
		
		//new method for calculating excitonic interactions
		//this uses the binding energy for excitons constrained to NPs, calculated under the effective mass approximation
		double C_coul = 1.786; //for approximately spherical nanoparticles
		double excitonBE_SP = (sourceNP.occupationTotalElectron <= sourceNP.occupationTotalHoles) ? 
				-C_coul*2*2/(sourceNP.dcin*sourceNP.diameter) : 0; //note: e^2 = 2 in Rydberg units
		double excitonBE_TP = (targetNP.occupationTotalElectron + 1 <= targetNP.occupationTotalHoles) ? 
				-C_coul*2*2/(targetNP.dcin*targetNP.diameter) : 0; //note: e^2 = 2 in Rydberg units
				
		return excitonBE_TP - excitonBE_SP;
    }

	private double nearestNeighborPoisson(Nanoparticle targetNP, Nanoparticle sourceNP, Sample sample) {
		double nnPoisson=0;
		double ccdistance; 
	    // for the electron before hopping this comes with negative sign as this is subtracted
		for(Nanoparticle sourceNeighbor: sourceNP.nearestNeighbors){
			ccdistance = sourceNP.centerDistanceMap.get(sourceNeighbor);
			nnPoisson += -Constants.e2*(sourceNeighbor.occupationTotalElectron-sourceNeighbor.occupationTotalHoles) / sample.dcout / ccdistance;
		}
		// for the electron after hopping
		for(Nanoparticle targetNeighbor: targetNP.nearestNeighbors){
			ccdistance = targetNP.centerDistanceMap.get(targetNeighbor);
			
			if(targetNeighbor==sourceNP)
				nnPoisson += Constants.e2*(targetNeighbor.occupationTotalElectron - targetNeighbor.occupationTotalHoles -1) / sample.dcout / ccdistance;	
			else
				nnPoisson += Constants.e2*(targetNeighbor.occupationTotalElectron- targetNeighbor.occupationTotalHoles) / sample.dcout / ccdistance;
		}
		return nnPoisson;
	}

    private double calculateRate(Nanoparticle targetNP, Nanoparticle sourceNP, int targetOrbital, int sourceOrbital, Sample sample){
    	double energy_diff = 0.0;
    	double rate = 0, overlap;
    	double npdistance = targetNP.edgeDistanceMap.get(sourceNP);
    	// Kinetic energy difference
    	energy_diff += targetNP.cbenergy[targetOrbital] - sourceNP.cbenergy[sourceOrbital] ;
    	
    	// On-Site charging energy
    	energy_diff += calculateCharging(targetNP, sourceNP, sample);
    	
    	// Electron-hole exciton interaction
    	energy_diff += calculateExciton(targetNP, sourceNP, sample);
    	//if(calculateExciton(targetNP, sourceNP, sample) > 0.0) System.out.println("exciton was used");
    	
    	// External potential
        energy_diff += -1.0*Constants.sqrt2*(sample.voltage/sample.cellz)*(targetNP.z-sourceNP.z-sample.cellz*Math.round((targetNP.z-sourceNP.z)/sample.cellz));

    	// Nearest neighbor interaction
        if(sample.poissonSolver=="nn")
        	energy_diff += nearestNeighborPoisson(targetNP, sourceNP, sample);
        
        overlap = Math.sqrt(-sample.emass*(sourceNP.cbenergy[sourceOrbital] + targetNP.cbenergy[targetOrbital]) / 2.0) ;

        // Two cases: Miller-Abrahms or Marcus
        switch (sample.hoppingMechanism) {
		case "ma":
			if (npdistance > 0) rate = sample.jumpfreq * Math.exp(-2.0*npdistance*overlap); else rate = sample.jumpfreq;
			if(energy_diff>0.0)
				rate = rate*Math.exp(-energy_diff/sample.temperature);
			break;

		case "marcus":
			rate = Constants.tpi*sample.marcusprefac2*Math.exp(-2.0 * npdistance * overlap) / Math.sqrt(Constants.fpi*sample.reorgenergy*sample.temperature);
	        rate = rate*Math.exp(-Math.pow((sample.reorgenergy + energy_diff), 2.0) / (4.0*sample.reorgenergy*sample.temperature));
			break;
		}
    	
    	return rate;
    }
    
    
    public void move(Nanoparticle sourceNP, int sourceBand, Nanoparticle destinationNP, int destinationBand) {
		
    	//System.out.println("moving electron "+this);
    	//System.out.println("destinationNP is "+destinationNP+", sourceNP is "+sourceNP);
		// remove from old NP first
    	sourceNP.remove_electron(this, sourceBand);
    	// add to the new NP
    	this.hostNP = destinationNP;
    	this.hostOrbital = destinationBand;
    	destinationNP.add_electron(this, destinationBand);
    	
    	double[] sourceDestination = {destinationNP.x, destinationNP.z};
    	this.visitedNPs.add(sourceDestination);
	}
    
    
    public double calculateEnergy() {
    	return 0;
		
	}


}
