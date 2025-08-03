#ifndef PhotonFlux_h
#define PhotonFlux_h


#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>

#include "gsl/gsl_sf_bessel.h"


using namespace std;

class PhotonFlux {

    public:
        PhotonFlux();
        ~PhotonFlux();

        double getDreesZeppenfeld(double x);
        double getModifiedDZ(double x);
        double getA(double x);
        double getKniehl(double x);
        double getJackson(double x);

        double xFromRapidity(double y){ return m_VM * exp(y)/squareRootOfS; }
        double xFromOmega(double omega){ return 2*omega/squareRootOfS; }

        double getDreesZeppenfeldFromRapidity(double y) { return getDreesZeppenfeld( xFromRapidity(y) ); }
        double getModifiedDZFromRapidity(double y) { return getModifiedDZ( xFromRapidity(y) ); }
        double getKniehlFromRapidity(double y) { return getKniehl( xFromRapidity(y) ); }
        double getJacksonFromRapidity(double y) { return getJackson( xFromRapidity(y) ); }

        double getDreesZeppenfeldFromPhotonEnergy(double omega) {double x = xFromOmega(omega); return x > 1 ? 0 : getDreesZeppenfeld(x); }
        double getModifiedDZFromPhotonEnergy(double omega) { double x = xFromOmega(omega); return x > 1 ? 0 : getModifiedDZ(x); }
        double getKniehlFromPhotonEnergy(double omega) { double x = xFromOmega(omega); return x > 1 ? 0 : getKniehl(x); }
        double getJacksonFromPhotonEnergy(double omega) { double x = xFromOmega(omega); return x > 1 ? 0 : getJackson(x); }
    private:
        double squareRootOfS = 510; //GeV
        double m_VM = 3.1;  //GeV
        double M_p = 0.939; //mass of proton in GeV  
        double alpha = 1.0/137; // fine structure constant

        vector<double> KniehlConstantsC = {-0.0276, 3.96, 13.8, -2.48, -0.891, -11.3, -0.716, -4.43, 0.238, -2.12};  

        double R_p = 0.77;       //fm
        double b_min = 2*R_p;    //fm
        int Z = 82;


};


#endif