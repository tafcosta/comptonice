/*
 * Photon.h
 *
 *  Created on: 12 Aug 2024
 *      Author: ntc132
 */

#ifndef SRC_PHOTON_H_
#define SRC_PHOTON_H_

#include "SimulationDependencies.h"

class Photon {
public:
    Photon(PhotonSpectrum &photonSpectrum);

	virtual ~Photon();

	double energy;
	double energyInitial;
	bool escaped;

	std::vector<double> direction;
	std::vector<double> position;

protected:
	PhotonSpectrum &photonSpectrum;

};

#endif /* SRC_PHOTON_H_ */
