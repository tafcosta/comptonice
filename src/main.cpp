#include "SimulationDependencies.h"

Grid *grid = new Grid(0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1, 256);
PhotonSpectrum *photonSpectrum = new PhotonSpectrumPowerLaw();
Scatter *scatter = new ScatterCompton();

static int nPhotons = 100000;

int main(){
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> RandomNumberZeroToOne(0.0, 1.0);

    std::ofstream outputFile("FinalOutput.txt");

    std::vector<Photon*> photons;
    photons.reserve(nPhotons);

    for(int iPhoton = 0; iPhoton < nPhotons; iPhoton++)
    	photons.emplace_back(new Photon(*grid, *photonSpectrum, iPhoton));

    std::cout << "Start Comptonisation" << std::endl;

    size_t totalPhotons = photons.size();
    size_t counter = 0;
    int nextProgress = 5; // Start with the first progress threshold at 5%

    for(Photon* photon : photons){

    	while(photon->insideDomain){
    		double opticalDepth = -std::log(1 - RandomNumberZeroToOne(gen));

    		photon->propagate(opticalDepth);

    		if(photon->insideDomain){
    			Electron *electron = new ElectronMaxwellBoltzmann(*grid);
    			scatter->doScattering(electron, photon);
    			delete electron;
    		}
    	}

    	counter++;
        float progress = (static_cast<float>(counter) / totalPhotons) * 100;
        if (progress >= nextProgress) {
            std::cout << "Progress: " << nextProgress << "%\r" << std::flush;
            nextProgress += 5; // Set the next progress threshold
        }

	}

    for(Photon* photon : photons){
        outputFile << photon->photonIndex << " " << photon->thetaLabFrame << " " << photon->energy << " " << photon->energyInitial << " " << photon->nScatters << std::endl;
    }

    outputFile.close();

    std::cout << "Comptonisation done." << std::endl;

    for(Photon* photon : photons)
    	delete photon;

    delete grid;
    delete photonSpectrum;
    delete scatter;
	return 0;
}
