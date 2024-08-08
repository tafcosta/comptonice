#include "SimulationDependencies.h"

double getElectronSpeed(double gasTemperatureKeV, double randomNumber);
double MaxwellBoltzmannCumulativeDistribution(double beta, double gasTemperatureKeV);
static double electronRestMassEnergyKeV = 511.;
static double MinBeta = 0.;
static double MaxBeta = 1.;
static double TOL = 0.0001;

int main(){

	double gasTemperatureKeV = 0.001;

	double photonPosition[3] = {0.,0.,0.};
	double photonDir[3] = {0.,0.,0.};
	double electronDir[3] = {0.,0.,0.};
	double electronSpeed = 0.;

	double photonEnergy = 10.0;
	int propagatePhoton = 1;
	int numberScatterings = 0;
	double randomNumber = 1;
	double criticalOpticalDepth, phi, theta, scatteringAngle;



    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> DistributionZeroToOne(0.0, 1.0);
    std::uniform_real_distribution<> azimuth(0.0, 2 * M_PI);
    std::uniform_real_distribution<> inclination(-1., 1.);


	while(propagatePhoton){

		phi   = azimuth(gen);
		theta = std::acos(inclination(gen));

		electronDir[0] = std::sin(theta) * std::cos(phi);
		electronDir[1] = std::sin(theta) * std::sin(phi);
		electronDir[2] = std::cos(theta);

		randomNumber  = DistributionZeroToOne(gen);
		electronSpeed = getElectronSpeed(gasTemperatureKeV, randomNumber);




		criticalOpticalDepth = -std::log(1 - DistributionZeroToOne(gen));
		phi   = azimuth(gen);
		theta = std::acos(inclination(gen));

		scatteringAngle = std::sin(theta) * std::cos(phi) * photonDir[0] + std::sin(theta) * std::sin(phi) * photonDir[1] +  std::cos(theta) * photonDir[2];
		photonEnergy = photonEnergy / (1 + photonEnergy/electronRestMassEnergyKeV * (1 - std::cos(scatteringAngle)));

		photonDir[0] = std::sin(theta) * std::cos(phi);
		photonDir[1] = std::sin(theta) * std::sin(phi);
		photonDir[2] = std::cos(theta);

		photonPosition[0] += criticalOpticalDepth * photonDir[0];
		photonPosition[1] += criticalOpticalDepth * photonDir[1];
		photonPosition[2] += criticalOpticalDepth * photonDir[2];

		numberScatterings += 1;

		if((photonPosition[0]*photonPosition[0] + photonPosition[1]*photonPosition[1] + photonPosition[2]*photonPosition[2]) > 100.)
			break;
	}

	std::cout << photonEnergy << " " << numberScatterings << std::endl;

	return 0;
}

double getElectronSpeed(double gasTemperatureKeV, double randomNumber){
	double electronSpeed;

	for(int nIter = 0; nIter < 100; nIter++){
		double lowerLimit = MaxwellBoltzmannCumulativeDistribution(MinBeta, gasTemperatureKeV) - randomNumber;
		double midPoint   = MaxwellBoltzmannCumulativeDistribution((MinBeta + MaxBeta)/2., gasTemperatureKeV) - randomNumber;

		if((midPoint == 0) || (MaxBeta - MinBeta)/2. < TOL){
			electronSpeed = (MinBeta + MaxBeta)/2.;
			break;
		}

		if(lowerLimit * midPoint > 0)
			MinBeta = (MinBeta + MaxBeta)/2.;
		else
			MaxBeta = (MinBeta + MaxBeta)/2.;
	}

	return electronSpeed;
}

double MaxwellBoltzmannCumulativeDistribution(double beta, double gasTemperatureKeV){
	return std::erf(electronRestMassEnergyKeV/gasTemperatureKeV/2. * beta * beta) - beta * std::sqrt(2/M_PI * electronRestMassEnergyKeV/gasTemperatureKeV) * std::exp(-electronRestMassEnergyKeV/gasTemperatureKeV/2. * beta * beta);
}
