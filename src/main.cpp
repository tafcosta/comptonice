#include "SimulationDependencies.h"

Grid *grid = new Grid(0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1, 128);
PhotonSpectrum *photonSpectrum = new PhotonSpectrumPowerLaw();
Scatter *scatter = new ScatterCompton();

static int nPhotons = 50000;

int main(){
	double gasTemperatureKeV = 0.001;
	double coupledEnergy = 0.;

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> RandomNumberZeroToOne(0.0, 1.0);

    std::ofstream outputFile("outputEnergy.txt");

    std::vector<Photon*> photons;
    photons.reserve(nPhotons);

    for(int iPhoton = 0; iPhoton < nPhotons; iPhoton++)
    	photons.emplace_back(new Photon(*grid, *photonSpectrum));

    std::cout << "Start Comptonisation" << std::endl;

    for(Photon* photon : photons){
    	while(photon->insideDomain){
    		double opticalDepth = -std::log(1 - RandomNumberZeroToOne(gen));

    		photon->propagate(opticalDepth);

    		if(!photon->insideDomain)
    			break;

    		scatter->doScattering(photon);
    	}
    	coupledEnergy += 1 - photon->energy/photon->energyInitial;
	}

    for(Photon* photon : photons)
        outputFile << photon->thetaLabFrame << " " << photon->phiLabFrame << std::endl;

    outputFile.close();

    std::cout << "Comptonisation done, energy transferred from photons = " << coupledEnergy/nPhotons << std::endl;

    for(Photon* photon : photons)
    	delete photon;

    delete grid;
    delete photonSpectrum;
    delete scatter;
	return 0;
}
