/*
 * PhotonSpectrumPowerLaw.h
 *
 *  Created on: 12 Aug 2024
 *      Author: ntc132
 */

#ifndef SRC_PHOTONSPECTRUMPOWERLAW_H_
#define SRC_PHOTONSPECTRUMPOWERLAW_H_

#include "SimulationDependencies.h"

class PhotonSpectrumPowerLaw: public PhotonSpectrum {
public:
	PhotonSpectrumPowerLaw();
	virtual ~PhotonSpectrumPowerLaw();

	double setPhotonEnergy() override;
};

#endif /* SRC_PHOTONSPECTRUMPOWERLAW_H_ */
