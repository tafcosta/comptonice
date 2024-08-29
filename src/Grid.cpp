/*
 * Grid.cpp
 *
 *  Created on: 20 Aug 2024
 *      Author: ntc132
 */

#include "SimulationDependencies.h"

Grid::Grid(double xmin, double xmax, double ymin, double ymax, double zmin, double zmax, int nGhost, int nCells) : xmin(xmin), ymin(ymin), zmin(zmin), nGhost(nGhost), nCells(nCells) {

	quantities =  std::vector<std::vector<std::vector<std::vector<int> > > > (nCells + 2*nGhost,
			std::vector<std::vector<std::vector<int> > >(nCells + 2*nGhost,
			std::vector<std::vector<int> >(nCells + 2*nGhost,
			std::vector<int>(nvar, 1.e-2))));

	dx = (xmax - xmin)/nCells;
	dy = (ymax - ymin)/nCells;
	dz = (zmax - zmin)/nCells;

	minXIndex = nGhost;
	maxXIndex = nGhost + nCells - 1;

	minYIndex = nGhost;
	maxYIndex = nGhost + nCells - 1;

	minZIndex = nGhost;
	maxZIndex = nGhost + nCells - 1;

	for(int i = minXIndex; i <= maxXIndex; i++){
		double x = getX(i) - 0.5;

		for(int j = minYIndex; j <= maxYIndex; j++){
			double y = getX(j) - 0.5;

			for(int k = minZIndex; k <= maxZIndex; k++){
				double z = getX(k) - 0.5;

				if((std::sqrt(x*x + y*y) < 0.4) && (std::fabs(z) < 0.05))
					quantities[i][j][k][DENS] = 2.e2;

			}
		}
	}
}

Grid::~Grid() {
	// TODO Auto-generated destructor stub
}
