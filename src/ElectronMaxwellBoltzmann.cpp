/*
 * ElectronMaxwellBoltzmann.cpp
 *
 *  Created on: 13 Aug 2024
 *      Author: ntc132
 */

#include "SimulationDependencies.h"

ElectronMaxwellBoltzmann::ElectronMaxwellBoltzmann(Grid &grid) : Electron(grid){

}

void ElectronMaxwellBoltzmann::getElectronSpeed(Photon* photon){
	double electronSpeed = 0.;
	double MinBeta  = 0.;
	double MaxBeta  = 1.;

    std::random_device rd;
    std::uniform_real_distribution<> DistributionZeroToOne(0.0, 1.0);
    std::mt19937 gen(rd());

    std::array<int, 3> indices = grid.getCellIndex(photon->position);
    double gasTemperatureKeV = grid.quantities[indices[0]][indices[1]][indices[2]][grid.TEMP];

    double randomNumber = DistributionZeroToOne(gen);

	for(int nIter = 0; nIter < nIterMax; nIter++){
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

	speed = electronSpeed;
	gamma = 1.0 / std::sqrt(1 - speed * speed);
}


double ElectronMaxwellBoltzmann::MaxwellBoltzmannCumulativeDistribution(double beta, double gasTemperatureKeV){
	return std::erf(electronRestMassEnergyKeV/gasTemperatureKeV/2. * beta * beta) - beta * std::sqrt(2/M_PI * electronRestMassEnergyKeV/gasTemperatureKeV) * std::exp(-electronRestMassEnergyKeV/gasTemperatureKeV/2. * beta * beta);
}

ElectronMaxwellBoltzmann::~ElectronMaxwellBoltzmann() {
	// TODO Auto-generated destructor stub
}

