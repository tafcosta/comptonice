/*
 * PhotonSpectrumPowerLaw.cpp
 *
 *  Created on: 12 Aug 2024
 *      Author: ntc132
 */

#include "SimulationDependencies.h"

void PhotonSpectrumPowerLaw::setInitialPhotonEnergy() {

		double minEnergy = 4.;
		double maxEnergy = 14.;

	    std::random_device rd;
	    std::mt19937 gen(rd());
	    std::uniform_real_distribution<> DistributionZeroToOne(0.0, 1.0);

	    for(int iPhoton = 0; iPhoton < photon.nPhotons; iPhoton++)
	    	photon.energyInitial[iPhoton] = minEnergy * std::pow(maxEnergy/minEnergy, DistributionZeroToOne(gen));
}


PhotonSpectrumPowerLaw::~PhotonSpectrumPowerLaw() {
	// TODO Auto-generated destructor stub
}

