/*
 * Electron.h
 *
 *  Created on: 13 Aug 2024
 *      Author: ntc132
 */

#ifndef SRC_ELECTRON_H_
#define SRC_ELECTRON_H_

class Electron {
public:
	Electron();
	virtual ~Electron();

	virtual double getElectronSpeed(double gasTemperatureKeV){};

	double electronRestMassEnergyKeV = 511.0;
};

#endif /* SRC_ELECTRON_H_ */
