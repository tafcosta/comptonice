/*
 * Electron.cpp
 *
 *  Created on: 13 Aug 2024
 *      Author: ntc132
 */

#include "SimulationDependencies.h"
#include "Grid.h"

Electron::Electron(Grid &grid) : grid(grid), direction(3, 0.0), speed(0.0), gamma(0.0) {

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> azimuth(0.0, 2 * M_PI);
    std::uniform_real_distribution<> cosInclination(-1., 1.);

	phiLabFrame = azimuth(gen);
	thetaLabFrame = cosInclination(gen);

	direction[0] = std::sin(thetaLabFrame) * std::cos(phiLabFrame);
	direction[1] = std::sin(thetaLabFrame) * std::sin(phiLabFrame);
	direction[2] = std::cos(thetaLabFrame);

}

Electron::~Electron() {

}

