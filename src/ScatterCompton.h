/*
 * ScatterCompton.h
 *
 *  Created on: 13 Aug 2024
 *      Author: ntc132
 */

#ifndef SRC_SCATTERCOMPTON_H_
#define SRC_SCATTERCOMPTON_H_

#include "SimulationDependencies.h"

class ScatterCompton: public Scatter {
public:
	virtual ~ScatterCompton();

	void doScattering(Photon*) override;

	double electronRestMassEnergyKeV = 511.0;

	void energyChangeCompton(Photon* photon, double scatteringAngle);
	void getPhotonDirectionLabFrame(double phiPhotonFrame, double thetaPhotonFrame, Photon* photon);
	double getScatteringAnglePhotonFrame(double photonEnergy);
	double KleinNishinaCumulativeDistribution(double photonEnergy, double theta);

protected:
	double TOL = 1.e-6;
	int nIterMax = 100;
};

#endif /* SRC_SCATTERCOMPTON_H_ */
