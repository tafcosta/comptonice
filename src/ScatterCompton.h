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

	double energyChangeCompton(double photonEnergy, double scatteringAngle);
	double getScatteringAngleLabFrame(double thetaLab, double scatteringAngle);
	double getScatteringAnglePhotonFrame(double photonEnergy);
	double KleinNishinaCumulativeDistribution(double photonEnergy, double theta);

protected:
	double TOL = 0.0001;
};

#endif /* SRC_SCATTERCOMPTON_H_ */
