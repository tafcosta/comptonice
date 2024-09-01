/*
 * Electron.h
 *
 *  Created on: 13 Aug 2024
 *      Author: ntc132
 */

#include "SimulationDependencies.h"
#include "Grid.h"

#ifndef SRC_ELECTRON_H_
#define SRC_ELECTRON_H_

class Electron {
public:
	Electron(Grid &grid);
	virtual ~Electron();
	virtual void getElectronSpeed(Photon* photon){};

	std::vector<double> direction;
	double phiLabFrame;
	double thetaLabFrame;
	double speed;
	double gamma;

	static const double electronRestMassEnergyKeV = 511.0;

protected:
	Grid &grid;
};

#endif /* SRC_ELECTRON_H_ */
