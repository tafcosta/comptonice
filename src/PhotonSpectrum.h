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
	PhotonSpectrum(Photon &photon);
	virtual ~PhotonSpectrum();

	virtual void setInitialPhotonEnergy(){};

protected:
	Photon &photon;
};

#endif /* SRC_PHOTONSPECTRUM_H_ */
