#include "PhotonFlux.h"
#include "matplotlibcpp.h"
#include <cmath>
#include <vector>

using namespace std;
namespace plt = matplotlibcpp;

void plotFractionOfMomentum();
void plotRapidity();
void plotPhotonEnergy();


int main(int argc, char* argv[]){

    plotFractionOfMomentum();
    plotRapidity();
    plotPhotonEnergy();


    return 0;
}

void plotRapidity(){
    // initialize PhotonFlux
    PhotonFlux* photonFlux = new PhotonFlux();

    vector<double> Y, DZ, DZM, K, J;
    int nPoints = 1000;
    double yMin = -4.0, yMax = 4.0;
    for(int i = 0; i < nPoints; ++i){
        double y = yMin + (yMax - yMin) * i / (nPoints - 1);
        Y.push_back(y);
        DZ.push_back(photonFlux->getDreesZeppenfeldFromRapidity(y));
        DZM.push_back(photonFlux->getModifiedDZFromRapidity(y));
        K.push_back(photonFlux->getKniehlFromRapidity(y));
        J.push_back(photonFlux->getJacksonFromRapidity(y));
    }

    plt::figure();
    plt::named_plot("Drees-Zeppenfeld",Y, DZ);
    plt::named_plot("Modified Drees-Zeppenfeld",Y, DZM);
    plt::named_plot("Kniehl",Y, K);
    plt::named_plot("Jackson",Y, J);
    PyRun_SimpleString("import matplotlib.pyplot as plt; plt.yscale('log')");
    plt::xlabel("y");
    plt::ylabel("dN/dy");
    plt::legend();
    plt::show();

    return;
}


void plotFractionOfMomentum(){
    // initialize PhotonFlux
    PhotonFlux* photonFlux = new PhotonFlux();

    vector<double> X, DZ, DZM, K, J;
    int nPoints = 1000;
    double xMin = 0.001, xMax = 0.2;
    for(int i = 0; i < nPoints; ++i){
        double x = xMin + (xMax - xMin) * i / (nPoints - 1);
        X.push_back(x);
        DZ.push_back(x*photonFlux->getDreesZeppenfeld(x));
        DZM.push_back(x*photonFlux->getModifiedDZ(x));
        K.push_back(x*photonFlux->getKniehl(x));
        J.push_back(x*photonFlux->getJackson(x));
    }

    plt::figure();
    plt::named_plot("Drees-Zeppenfeld",X, DZ);
    plt::named_plot("Modified Drees-Zeppenfeld",X, DZM);
    plt::named_plot("Kniehl",X, K);
    plt::named_plot("Jackson",X, J);
    PyRun_SimpleString("import matplotlib.pyplot as plt; plt.yscale('log')");
    plt::xlabel("x");
    plt::ylabel("xf(x)");
    plt::legend();
    plt::show();

    return;
}

void plotPhotonEnergy(){
    // initialize PhotonFlux
    PhotonFlux* photonFlux = new PhotonFlux();

    vector<double> Omega, DZ, DZM, K, J;
    int nPoints = 1000;
    double omegaMin = 0.01, omegaMax = 250;  //GeV
    for(int i = 0; i < nPoints; ++i){
        double omega = omegaMin + (omegaMax - omegaMin) * i / (nPoints - 1);
        Omega.push_back(omega);
        DZ.push_back(photonFlux->getDreesZeppenfeldFromPhotonEnergy(omega));
        DZM.push_back(photonFlux->getModifiedDZFromPhotonEnergy(omega));
        K.push_back(photonFlux->getKniehlFromPhotonEnergy(omega));
        J.push_back(photonFlux->getJacksonFromPhotonEnergy(omega));
    }

    plt::figure();
    plt::named_plot("Drees-Zeppenfeld",Omega, DZ);
    plt::named_plot("Modified Drees-Zeppenfeld",Omega, DZM);
    plt::named_plot("Kniehl",Omega, K);
    plt::named_plot("Jackson",Omega, J);
    PyRun_SimpleString("import matplotlib.pyplot as plt; plt.yscale('log')");
    plt::xlabel("ω (GeV)");
    plt::ylabel("dN/dω");
    plt::legend();
    plt::show();

    return;
}