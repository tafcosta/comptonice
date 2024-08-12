/*
 * PhotonSpectrum.h
 *
 *  Created on: 12 Aug 2024
 *      Author: ntc132
 */

#ifndef SRC_PHOTONSPECTRUM_H_
#define SRC_PHOTONSPECTRUM_H_

class PhotonSpectrum {
public:
	PhotonSpectrum();
	virtual ~PhotonSpectrum();

	virtual double setInitialPhotonEnergy(){return 0.;};
};

#endif /* SRC_PHOTONSPECTRUM_H_ */
