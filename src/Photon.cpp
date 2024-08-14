/*
 * Photon.cpp
 *
 *  Created on: 12 Aug 2024
 *      Author: ntc132
 */

#include "SimulationDependencies.h"


Photon::Photon(int nPhotons) : nPhotons(nPhotons) {

	energy = std::vector<double> (nPhotons, 0.0) ;
	energyInitial = std::vector<double> (nPhotons, 0.0) ;

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> azimuth(0.0, 2 * M_PI);
    std::uniform_real_distribution<> inclination(-1., 1.);

	for(int iPhoton = 0; iPhoton < nPhotons; iPhoton++){

		double phi      = azimuth(gen);
		double thetaLab = std::acos(inclination(gen));

		direction[iPhoton][0] = std::sin(thetaLab) * std::cos(phi);
		direction[iPhoton][1] = std::sin(thetaLab) * std::sin(phi);
		direction[iPhoton][2] = std::cos(thetaLab);
	}

	position = std::vector<std::vector<double> > (nPhotons, std::vector<double>(3, 0.0));

}

Photon::~Photon() {
	// TODO Auto-generated destructor stub
}

