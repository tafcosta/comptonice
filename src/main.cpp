#include "SimulationDependencies.h"

double energyChangeCompton(double photonEnergy, double scatteringAngle);
double getElectronSpeed(double gasTemperatureKeV);
double getScatteringAngleLabFrame(double thetaLab, double scatteringAngle);
double getScatteringAnglePhotonFrame(double photonEnergy);
double KleinNishinaCumulativeDistribution(double photonEnergy, double theta);
double MaxwellBoltzmannCumulativeDistribution(double beta, double gasTemperatureKeV);

static double electronRestMassEnergyKeV = 511.0;
static double TOL = 0.0001;

Photon *photon = new Photon(30000);
PhotonSpectrum *photonSpectrum = new PhotonSpectrumPowerLaw(*photon);

int main(){

	double domainRadius = 10.;
	double gasTemperatureKeV = 0.001;
	double coupledEnergy = 0.;

	double electronDir[3] = {0.,0.,0.};
	std::vector<double> energyOut;
	double criticalOpticalDepth, phi, thetaLab;

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> RandomNumberZeroToOne(0.0, 1.0);
    std::uniform_real_distribution<> azimuth(0.0, 2 * M_PI);
    std::uniform_real_distribution<> inclination(-1., 1.);

    std::ofstream outputFile("outputEnergy.txt");

	photonSpectrum->setInitialPhotonEnergy();

    for(int iPhoton = 0; iPhoton < photon->nPhotons; iPhoton++){

    	photon->energy[iPhoton] = photon->energyInitial[iPhoton];

    	while(int propagatePhoton = 1){

    		criticalOpticalDepth = -std::log(1 - RandomNumberZeroToOne(gen));

    		/*
    		phi   = azimuth(gen);
    		theta = std::acos(inclination(gen));

    		electronDir[0] = std::sin(theta) * std::cos(phi);
    		electronDir[1] = std::sin(theta) * std::sin(phi);
    		electronDir[2] = std::cos(theta);

    		double electronSpeed = getElectronSpeed(gasTemperatureKeV);
    		 */

    		phi   = azimuth(gen);
    		double scatteringAngle = getScatteringAnglePhotonFrame(photon->energy[iPhoton]);
    		thetaLab = getScatteringAngleLabFrame(thetaLab, scatteringAngle);

    		photon->direction[iPhoton][0] = std::sin(thetaLab) * std::cos(phi);
    		photon->direction[iPhoton][1] = std::sin(thetaLab) * std::sin(phi);
    		photon->direction[iPhoton][2] = std::cos(thetaLab);

    		for(int i = 0; i < 3; i++)
    			photon->position[iPhoton][i] += criticalOpticalDepth * photon->direction[iPhoton][i];

    		if((photon->position[iPhoton][0]*photon->position[iPhoton][0] + photon->position[iPhoton][1]*photon->position[iPhoton][1] + photon->position[iPhoton][2]*photon->position[iPhoton][2]) > std::pow(domainRadius,2)){
        		energyOut.push_back(photon->energy[iPhoton]);
    			coupledEnergy += 1 - photon->energy[iPhoton]/photon->energyInitial[iPhoton];
    			break;
    		}

    		photon->energy[iPhoton] = energyChangeCompton(photon->energy[iPhoton], scatteringAngle);

    	}
	}

    for (const double& element : energyOut) {
        outputFile << element << " ";
    }

    outputFile.close();

    std::cout << "coupled energy = " << coupledEnergy/photon->nPhotons << std::endl;

    delete photon;
    delete photonSpectrum;
	return 0;
}

double energyChangeCompton(double photonEnergy, double scatteringAngle){
	return photonEnergy / (1 + photonEnergy/electronRestMassEnergyKeV * (1 - std::cos(scatteringAngle)));
}

double getElectronSpeed(double gasTemperatureKeV){
	double electronSpeed = 0.;
	double MinBeta  = 0.;
	double MaxBeta  = 1.;

    std::random_device rd;
    std::uniform_real_distribution<> DistributionZeroToOne(0.0, 1.0);
    std::mt19937 gen(rd());

    double randomNumber = DistributionZeroToOne(gen);

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

double getScatteringAngleLabFrame(double thetaLab, double scatteringAngle){
	return thetaLab + scatteringAngle;
}

double getScatteringAnglePhotonFrame(double photonEnergy){
	double ScatteringAngle = 0.;
	double MinAngle = 0.;
	double MaxAngle = M_PI;

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> DistributionZeroToOne(0.0, 1.0);

    double randomNumber = DistributionZeroToOne(gen);

	for(int nIter = 0; nIter < 100; nIter++){
		double lowerLimit = KleinNishinaCumulativeDistribution(photonEnergy, MinAngle) - randomNumber;
		double midPoint   = KleinNishinaCumulativeDistribution(photonEnergy, (MinAngle + MaxAngle)/2.) - randomNumber;

		if((midPoint == 0) || (MaxAngle - MinAngle)/2. < TOL){
			ScatteringAngle = (MinAngle+ MaxAngle)/2.;
			break;
		}

		if(lowerLimit * midPoint > 0)
			MinAngle = (MinAngle + MaxAngle)/2.;
		else
			MaxAngle = (MinAngle + MaxAngle)/2.;
	}

	return ScatteringAngle;

}

double KleinNishinaCumulativeDistribution(double photonEnergy, double theta){
	double x = photonEnergy/electronRestMassEnergyKeV;
	double norm = 3./4 * ((1 + x)/std::pow(x,3) * (2*x * (1+x)/(1 + 2*x) - std::log(1 + 2*x)) + std::log(1 + 2*x) /(2*x) - (1 + 3*x)/std::pow(1 + 2*x,2));

	return ((2 + 6*x + std::pow(x,2)) / (2. * std::pow(x,3)) \
			+ (-2*x * std::cos(theta) + (-2. - 6*x - 5 * std::pow(x,2) + 2*x * (1 + 2*x) * std::cos(theta)) / std::pow(1 + x - x*std::cos(theta),2) \
			+ 2 * (-2. - 2*x + std::pow(x,2)) * std::log(1 + x - x * std::cos(theta)))/(2. * std::pow(x,3))) * 3/8./norm;

}

double MaxwellBoltzmannCumulativeDistribution(double beta, double gasTemperatureKeV){
	return std::erf(electronRestMassEnergyKeV/gasTemperatureKeV/2. * beta * beta) - beta * std::sqrt(2/M_PI * electronRestMassEnergyKeV/gasTemperatureKeV) * std::exp(-electronRestMassEnergyKeV/gasTemperatureKeV/2. * beta * beta);
}
