/*
 * PhotonSpectrumSingleFrequency.cpp
 *
 *  Created on: 13 Aug 2024
 *      Author: ntc132
 */

#include "SimulationDependencies.h"

PhotonSpectrumSingleFrequency::PhotonSpectrumSingleFrequency(double photonEnergy): photonEnergy(photonEnergy) {

}

double PhotonSpectrumSingleFrequency::setPhotonEnergy(){
    return photonEnergy;
}

PhotonSpectrumSingleFrequency::~PhotonSpectrumSingleFrequency() {
	// TODO Auto-generated destructor stub
}

