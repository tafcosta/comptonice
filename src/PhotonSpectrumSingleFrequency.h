/*
 * PhotonSpectrumSingleFrequency.h
 *
 *  Created on: 13 Aug 2024
 *      Author: ntc132
 */

#ifndef SRC_PHOTONSPECTRUMSINGLEFREQUENCY_H_
#define SRC_PHOTONSPECTRUMSINGLEFREQUENCY_H_

#include "SimulationDependencies.h"

class PhotonSpectrumSingleFrequency: public PhotonSpectrum {
public:
	PhotonSpectrumSingleFrequency(double photonEnergy);
	virtual ~PhotonSpectrumSingleFrequency();

	double setPhotonEnergy() override;

	double photonEnergy;
};

#endif /* SRC_PHOTONSPECTRUMSINGLEFREQUENCY_H_ */
