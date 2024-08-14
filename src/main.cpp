#include "SimulationDependencies.h"

PhotonSpectrum *photonSpectrum = new PhotonSpectrumPowerLaw();
Scatter *scatter = new ScatterCompton();

int photonEscapes(double criticalOpticalDepth, std::vector<double> position);

static double domainRadius = 30.;
static int nPhotons = 30000;

int main(){
	double gasTemperatureKeV = 0.001;
	double coupledEnergy = 0.;

	std::vector<double> energyOut;
	double criticalOpticalDepth;

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> RandomNumberZeroToOne(0.0, 1.0);

    std::ofstream outputFile("outputEnergy.txt");

    std::vector<Photon*> photons;
    photons.reserve(nPhotons);

    std::cout << "Start Comptonisation" << std::endl;

    for(int iPhoton = 0; iPhoton < nPhotons; iPhoton++)
    	photons.emplace_back(new Photon(*photonSpectrum));

    for(Photon* photon : photons){

    	while(true){
    		criticalOpticalDepth = -std::log(1 - RandomNumberZeroToOne(gen));

    		if(photonEscapes(criticalOpticalDepth, photon->position)){
    			energyOut.push_back(photon->energy);
    			coupledEnergy += 1 - photon->energy/photon->energyInitial;

    			photon->escaped = true;
    			break;
    		}

    		scatter->doScattering(photon);

    		for(int i = 0; i < 3; i++)
    			photon->position[i] += criticalOpticalDepth * photon->direction[i];
    	}
	}

    for (const double& element : energyOut) {
        outputFile << element << " ";
    }

    outputFile.close();

    std::cout << "Comptonisation done, coupled energy = " << coupledEnergy/nPhotons << std::endl;

    for(Photon* photon : photons)
    	delete photon;

    delete photonSpectrum;
    delete scatter;
	return 0;
}

int photonEscapes(double criticalOpticalDepth, std::vector<double> position){
	if((position[0]*position[0] + position[1]*position[1] + position[2]*position[2]) > std::pow(domainRadius,2))
		return 1;

	return 0;
}
