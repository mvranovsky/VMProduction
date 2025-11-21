
#include "VMProduction.h"


int main(){

    cout << "Starting VMProduction plots..." << endl;
    cout << "C++ standard version: " << __cplusplus << endl;

    cout << "Generating plot of fraction of momentum..." << endl;
    plotFractionOfMomentum();
    cout << "Generating plot of rapidity..." << endl;
    plotRapidity();
    
    cout << "Generating plot of photon energy..." << endl;
    plotPhotonEnergy();
    
    cout << "Generating plot of dipole cross section..." << endl;
    plotDipoleCrossSection(0.0001, 0.001, 0.1, true);

    cout << "Finished generating all plots." << endl;
    return 0;
}
