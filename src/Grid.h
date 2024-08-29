/*
 * Grid.h
 *
 *  Created on: 20 Aug 2024
 *      Author: ntc132
 */

#include "SimulationDependencies.h"

#ifndef SRC_GRID_H_
#define SRC_GRID_H_

class Grid {
public:
	Grid(double xmin, double xmax, double ymin, double ymax, double zmin, double zmax, int nGhost, int nCells);
	virtual ~Grid();

	double dx, dy, dz;
	double xmin;
	double ymin;
	double zmin;

	int minXIndex, maxXIndex;
	int minYIndex, maxYIndex;
	int minZIndex, maxZIndex;

	int nGhost, nCells;

	int nvar = 2;

	std::vector<std::vector<std::vector<std::vector<int> > > > quantities;

	static const int DENS   = 0;
	static const int TEMP   = 1;

	std::array<int, 3> getCellIndex(std::vector<double> position){
		std::array<int, 3> indices;

		indices[0] = minXIndex + std::floor((position[0] - xmin)/dx);
		indices[1] = minYIndex + std::floor((position[1] - ymin)/dy);
		indices[2] = minZIndex + std::floor((position[2] - zmin)/dz);

		return indices;
	}

	double getX(int cellIndex){
		return xmin + (cellIndex - minXIndex) * dx + dx/2.;
	}
};

#endif /* SRC_GRID_H_ */

