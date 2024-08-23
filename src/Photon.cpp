/*
 * Photon.cpp
 *
 *  Created on: 12 Aug 2024
 *      Author: ntc132
 */

#include "SimulationDependencies.h"

Photon::Photon(Grid &grid, PhotonSpectrum &photonSpectrum) : grid(grid), photonSpectrum(photonSpectrum), insideDomain(true), direction(3, 0.0), position(3, 0.50001){

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

void Photon::propagate(double opticalDepth){
	std::array<double, 3> perpendicularVector;
	std::array<double, 3> pointOnPlane;
	bool propagate = true;

	while(propagate){
		std::array<int, 3> indices = grid.getCellIndex(position);
		double distanceToEdge = 1.e10;
		double remainingPhotonPath = opticalDepth/grid.quantities[indices[0]][indices[1]][indices[2]][grid.DENS];


		for(int dir = 0; dir < 6; dir++){

			if(dir == 0){
				perpendicularVector = perpendicularVectorXpositive;
				pointOnPlane[0] = 0.5 * (grid.getX(indices[0] + 1) + grid.getX(indices[0]));
				pointOnPlane[1] = grid.getX(indices[1]);
				pointOnPlane[2] = grid.getX(indices[2]);
			}

			if(dir == 1){
				perpendicularVector = perpendicularVectorXnegative;
				pointOnPlane[0] = 0.5 * (grid.getX(indices[0] - 1) + grid.getX(indices[0]));
				pointOnPlane[1] = grid.getX(indices[1]);
				pointOnPlane[2] = grid.getX(indices[2]);
			}

			if(dir == 2){
				perpendicularVector = perpendicularVectorYpositive;
				pointOnPlane[0] = grid.getX(indices[0]);
				pointOnPlane[1] = 0.5 * (grid.getX(indices[1] + 1) + grid.getX(indices[1]));
				pointOnPlane[2] = grid.getX(indices[2]);
			}

			if(dir == 3){
				perpendicularVector = perpendicularVectorYnegative;
				pointOnPlane[0] = grid.getX(indices[0]);
				pointOnPlane[1] = 0.5 * (grid.getX(indices[1] - 1) + grid.getX(indices[1]));
				pointOnPlane[2] = grid.getX(indices[2]);
			}

			if(dir == 4){
				perpendicularVector = perpendicularVectorZpositive;
				pointOnPlane[0] = grid.getX(indices[0]);
				pointOnPlane[1] = grid.getX(indices[1]);
				pointOnPlane[2] = 0.5 * (grid.getX(indices[2] + 1) + grid.getX(indices[2]));
			}

			if(dir == 5){
				perpendicularVector = perpendicularVectorZnegative;
				pointOnPlane[0] = grid.getX(indices[0]);
				pointOnPlane[1] = grid.getX(indices[1]);
				pointOnPlane[2] = 0.5 * (grid.getX(indices[2] - 1) + grid.getX(indices[2]));
			}

			double distanceTmp  = perpendicularVector[0] * (pointOnPlane[0] - position[0]) + perpendicularVector[1] * (pointOnPlane[1] - position[1]) + perpendicularVector[2] * (pointOnPlane[2] - position[2]);
			distanceTmp /= (perpendicularVector[0] * direction[0] + perpendicularVector[1] * direction[1] + perpendicularVector[2] * direction[2]);

			if((distanceTmp > 0) && (distanceTmp < distanceToEdge))
				distanceToEdge = distanceTmp;
		}

		if(remainingPhotonPath < distanceToEdge){
			for(int k = 0; k < 3; k++)
				position[k] += direction[k] * remainingPhotonPath;

			propagate = false;
		} else{
			for(int k = 0; k < 3; k++)
				position[k] += direction[k] * distanceToEdge;

			opticalDepth -= distanceToEdge * grid.quantities[indices[0]][indices[1]][indices[2]][grid.DENS];
		}

		indices = grid.getCellIndex(position);
		if((indices[0] >= grid.maxXIndex) || (indices[0] <= grid.minXIndex) || (indices[1] >= grid.maxYIndex) || (indices[1] <= grid.minYIndex) || (indices[2] >= grid.maxZIndex) || (indices[2] <= grid.minZIndex)){
			propagate = false;
			insideDomain = false;
		}

	}


}

Photon::~Photon() {
	// TODO Auto-generated destructor stub
}

