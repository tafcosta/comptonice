/*
 * ElectronMaxwellBoltzmann.h
 *
 *  Created on: 13 Aug 2024
 *      Author: ntc132
 */

#ifndef SRC_ELECTRONMAXWELLBOLTZMANN_H_
#define SRC_ELECTRONMAXWELLBOLTZMANN_H_

#include "SimulationDependencies.h"
#include "Electron.h"
#include "Grid.h"

class ElectronMaxwellBoltzmann: public Electron {
public:
	ElectronMaxwellBoltzmann(Grid &grid);
	virtual ~ElectronMaxwellBoltzmann();

	void getElectronSpeed(Photon* photon) override;
	double MaxwellBoltzmannCumulativeDistribution(double beta, double gasTemperatureKeV);

protected:
	double TOL = 1.e-5;
	int nIterMax = 100;
};

#endif /* SRC_ELECTRONMAXWELLBOLTZMANN_H_ */
