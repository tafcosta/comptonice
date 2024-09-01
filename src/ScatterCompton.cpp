/*
 * ScatterCompton.cpp
 *
 *  Created on: 13 Aug 2024
 *      Author: ntc132
 */

#include "SimulationDependencies.h"

void ScatterCompton::doScattering(Electron* electron, Photon* photon) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> azimuth(0.0, 2 * M_PI);

	electron->getElectronSpeed(photon);

	double cosBetweenPhotonElectron = electron->direction[0] * photon->direction[0] + electron->direction[1] * photon->direction[1] + electron->direction[2] * photon->direction[2];

	double phiPhotonFrame = azimuth(gen);
	double thetaPhotonFrame = getScatteringAnglePhotonFrame(electron->gamma * (1 + cosBetweenPhotonElectron * electron->speed) * photon->energy);

	getPhotonDirectionLabFrame(phiPhotonFrame, thetaPhotonFrame, photon);
	energyChangeCompton(photon, thetaPhotonFrame);
	getEnergyInLabeFrame(electron, photon, cosBetweenPhotonElectron);

	photon->nScatters += 1;
}

void ScatterCompton::energyChangeCompton(Photon* photon, double scatteringAngle){
	 photon->energy = photon->energy / (1 + photon->energy/electronRestMassEnergyKeV * (1 - std::cos(scatteringAngle)));
}

void ScatterCompton::getEnergyInLabeFrame(Electron* electron, Photon* photon, double cosBetweenPhotonElectron){
	photon->energy = photon->energy * electron->gamma * (1 - cosBetweenPhotonElectron * electron->speed);
}

void ScatterCompton::getPhotonDirectionLabFrame(double phiPhotonFrame, double thetaPhotonFrame, Photon* photon){
	photon->direction[0] = std::sin(thetaPhotonFrame) * (std::cos(photon->phiLabFrame) * std::cos(photon->thetaLabFrame) * std::cos(phiPhotonFrame) - std::sin(photon->phiLabFrame)*std::sin(phiPhotonFrame)) + std::cos(photon->phiLabFrame)*std::sin(photon->thetaLabFrame)*std::cos(thetaPhotonFrame);
	photon->direction[1] = std::sin(thetaPhotonFrame) * (std::sin(photon->phiLabFrame) * std::cos(photon->thetaLabFrame) * std::cos(phiPhotonFrame) + std::cos(photon->phiLabFrame)*std::sin(phiPhotonFrame)) + std::sin(photon->phiLabFrame)*std::sin(photon->thetaLabFrame)*std::cos(thetaPhotonFrame);
	photon->direction[2] = -std::sin(photon->thetaLabFrame) * std::cos(phiPhotonFrame) * std::sin(thetaPhotonFrame) + std::cos(photon->thetaLabFrame) * std::cos(thetaPhotonFrame);

	double radius = std::sqrt(std::pow(photon->direction[0],2) + std::pow(photon->direction[1],2));
	photon->thetaLabFrame = std::atan2(radius, photon->direction[2]);
	photon->phiLabFrame = std::atan2(photon->direction[1], photon->direction[0]);
}

double ScatterCompton::getScatteringAnglePhotonFrame(double photonEnergy){
	double ScatteringAngle = 0.;
	double MinAngle = 0.;
	double MaxAngle = M_PI;

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> DistributionZeroToOne(0.0, 1.0);

    double randomNumber = DistributionZeroToOne(gen);

	for(int nIter = 0; nIter < nIterMax; nIter++){
		double lowerLimit = KleinNishinaCumulativeDistribution(photonEnergy, MinAngle) - randomNumber;
		double midPoint   = KleinNishinaCumulativeDistribution(photonEnergy, (MinAngle + MaxAngle)/2.) - randomNumber;

		if((midPoint == 0) || (MaxAngle - MinAngle)/2. < TOL){
			ScatteringAngle = (MinAngle + MaxAngle)/2.;
			break;
		}

		if(lowerLimit * midPoint > 0)
			MinAngle = (MinAngle + MaxAngle)/2.;
		else
			MaxAngle = (MinAngle + MaxAngle)/2.;
	}

	return ScatteringAngle;

}

double ScatterCompton::KleinNishinaCumulativeDistribution(double photonEnergy, double theta){
	double x = photonEnergy/electronRestMassEnergyKeV;
	double norm = 3./4 * ((1 + x)/std::pow(x,3) * (2*x * (1+x)/(1 + 2*x) - std::log(1 + 2*x)) + std::log(1 + 2*x) /(2*x) - (1 + 3*x)/std::pow(1 + 2*x,2));

	return ((2 + 6*x + std::pow(x,2)) / (2. * std::pow(x,3)) \
			+ (-2*x * std::cos(theta) + (-2. - 6*x - 5 * std::pow(x,2) + 2*x * (1 + 2*x) * std::cos(theta)) / std::pow(1 + x - x*std::cos(theta),2) \
			+ 2 * (-2. - 2*x + std::pow(x,2)) * std::log(1 + x - x * std::cos(theta)))/(2. * std::pow(x,3))) * 3/8./norm;
}

ScatterCompton::~ScatterCompton() {
	// TODO Auto-generated destructor stub
}

