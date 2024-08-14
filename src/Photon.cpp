/*
 * Photon.cpp
 *
 *  Created on: 12 Aug 2024
 *      Author: ntc132
 */

#include "SimulationDependencies.h"

Photon::Photon(PhotonSpectrum &photonSpectrum) : photonSpectrum(photonSpectrum), escaped(false), direction(3, 0.0), position(3, 0.0){

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> azimuth(0.0, 2 * M_PI);
    std::uniform_real_distribution<> inclination(-1., 1.);

    double phi      = azimuth(gen);
	double thetaLab = std::acos(inclination(gen));

	energyInitial = photonSpectrum.setPhotonEnergy();
	energy = energyInitial;

	direction[0] = std::sin(thetaLab) * std::cos(phi);
	direction[1] = std::sin(thetaLab) * std::sin(phi);
	direction[2] = std::cos(thetaLab);
}

Photon::~Photon() {
	// TODO Auto-generated destructor stub
}

