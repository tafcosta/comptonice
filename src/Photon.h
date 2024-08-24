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
    Photon(Grid &grid, PhotonSpectrum &photonSpectrum);

	virtual ~Photon();

	double energy;
	double energyInitial;
	double thetaLabFrame;
	double phiLabFrame;
	bool insideDomain;

	std::vector<double> direction;
	std::vector<double> position;

	void propagate(double opticalDepth);

protected:
	Grid &grid;
	PhotonSpectrum &photonSpectrum;

private:
	std::array<double, 3> perpendicularVectorXpositive = {1.,0.,0.};
	std::array<double, 3> perpendicularVectorXnegative = {-1.,0.,0.};
	std::array<double, 3> perpendicularVectorYpositive = {0.,1.,0.};
	std::array<double, 3> perpendicularVectorYnegative = {0.,-1.,0.};
	std::array<double, 3> perpendicularVectorZpositive = {0.,0.,1.};
	std::array<double, 3> perpendicularVectorZnegative = {0.,0.,-1.};
};

#endif /* SRC_PHOTON_H_ */
