/*
 * Scatter.h
 *
 *  Created on: 13 Aug 2024
 *      Author: ntc132
 */

#ifndef SRC_SCATTER_H_
#define SRC_SCATTER_H_

#include "SimulationDependencies.h"

class Scatter {
public:
	Scatter();
	virtual ~Scatter();

	virtual void doScattering(Electron*, Photon*){};

};

#endif /* SRC_SCATTER_H_ */
