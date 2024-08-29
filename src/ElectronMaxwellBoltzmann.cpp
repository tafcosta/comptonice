/*
 * ElectronMaxwellBoltzmann.cpp
 *
 *  Created on: 13 Aug 2024
 *      Author: ntc132
 */

#include "SimulationDependencies.h"

ElectronMaxwellBoltzmann::ElectronMaxwellBoltzmann() {

}

double ElectronMaxwellBoltzmann::getElectronSpeed(double gasTemperatureKeV){
	double electronSpeed = 0.;
	double MinBeta  = 0.;
	double MaxBeta  = 1.;

    std::random_device rd;
    std::uniform_real_distribution<> DistributionZeroToOne(0.0, 1.0);
    std::mt19937 gen(rd());

    double randomNumber = DistributionZeroToOne(gen);

	for(int nIter = 0; nIter < 100; nIter++){
		double lowerLimit = MaxwellBoltzmannCumulativeDistribution(MinBeta, gasTemperatureKeV) - randomNumber;
		double midPoint   = MaxwellBoltzmannCumulativeDistribution((MinBeta + MaxBeta)/2., gasTemperatureKeV) - randomNumber;

		if((midPoint == 0) || (MaxBeta - MinBeta)/2. < TOL){
			electronSpeed = (MinBeta + MaxBeta)/2.;
			break;
		}

		if(lowerLimit * midPoint > 0)
			MinBeta = (MinBeta + MaxBeta)/2.;
		else
			MaxBeta = (MinBeta + MaxBeta)/2.;
	}

	return electronSpeed;
}


double ElectronMaxwellBoltzmann::MaxwellBoltzmannCumulativeDistribution(double beta, double gasTemperatureKeV){
	return std::erf(electronRestMassEnergyKeV/gasTemperatureKeV/2. * beta * beta) - beta * std::sqrt(2/M_PI * electronRestMassEnergyKeV/gasTemperatureKeV) * std::exp(-electronRestMassEnergyKeV/gasTemperatureKeV/2. * beta * beta);
}

ElectronMaxwellBoltzmann::~ElectronMaxwellBoltzmann() {
	// TODO Auto-generated destructor stub
}

