/*
 * PhotonSpectrumPowerLaw.cpp
 *
 *  Created on: 12 Aug 2024
 *      Author: ntc132
 */

#include "SimulationDependencies.h"

PhotonSpectrumPowerLaw::PhotonSpectrumPowerLaw(){
}

double PhotonSpectrumPowerLaw::setPhotonEnergy(){
	double minEnergy = 4.;
	double maxEnergy = 14.;

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> DistributionZeroToOne(0.0, 1.0);

    return minEnergy * std::pow(maxEnergy/minEnergy, DistributionZeroToOne(gen));
}

PhotonSpectrumPowerLaw::~PhotonSpectrumPowerLaw() {
	// TODO Auto-generated destructor stub
}

