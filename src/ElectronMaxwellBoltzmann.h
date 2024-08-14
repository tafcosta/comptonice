/*
 * ElectronMaxwellBoltzmann.h
 *
 *  Created on: 13 Aug 2024
 *      Author: ntc132
 */

#ifndef SRC_ELECTRONMAXWELLBOLTZMANN_H_
#define SRC_ELECTRONMAXWELLBOLTZMANN_H_

#include "Electron.h"

class ElectronMaxwellBoltzmann: public Electron {
public:
	ElectronMaxwellBoltzmann();
	virtual ~ElectronMaxwellBoltzmann();

	double getElectronSpeed(double gasTemperatureKeV) override;
	double MaxwellBoltzmannCumulativeDistribution(double beta, double gasTemperatureKeV);

protected:
	double TOL = 0.0001;
};

#endif /* SRC_ELECTRONMAXWELLBOLTZMANN_H_ */
