/*
 * PhotonSpectrum.h
 *
 *  Created on: 12 Aug 2024
 *      Author: ntc132
 */

#ifndef SRC_PHOTONSPECTRUM_H_
#define SRC_PHOTONSPECTRUM_H_

#include "SimulationDependencies.h"

class PhotonSpectrum {
public:
	PhotonSpectrum();
	virtual ~PhotonSpectrum();

	double virtual setPhotonEnergy(){};
};

#endif /* SRC_PHOTONSPECTRUM_H_ */
