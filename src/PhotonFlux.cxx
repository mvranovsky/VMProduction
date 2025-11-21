#include <PhotonFlux.h>



PhotonFlux::PhotonFlux() {

}

PhotonFlux::~PhotonFlux() {
    // Destructor implementation if needed
    // Clean up resources if any were allocated
}

double PhotonFlux::getA(double x) {
    // factor A which comes up in Drees, Zeppenfeld Photon flux

    double Q2min = pow(x*M_p, 2)/(1-x);   //approximation used in https://doi.org/10.48550/arXiv.hep-ph/0412096
    return (1 + (0.71)/Q2min);
}



double PhotonFlux::getDreesZeppenfeld(double x){
    // Drees-Zeppenfeld:
    // https://doi.org/10.1103/PhysRevD.39.2536

    double A = getA(x);
    double val = alpha/(2*M_PI*x) *( 1 + pow(1-x, 2) )*(log(A) - 11/6 + 3/A - 3/( 2 * pow(A,2) ) + 1/(3*pow( A, 3)));
    return val;

}

double PhotonFlux::getModifiedDZ(double x){
    // modified Drees-Zeppenfeld, as described in http://arxiv.org/abs/hep-ph/0412096v1
    // when no approximation with term (Q^2 - Q^{2}_{min})/Q^4 #approx Q^{-2} is used

    double A = getA(x);
    double val = alpha/(2*M_PI*x) *( 1 + pow(1-x, 2) )*((A + 3)*log(A)/(A - 1) - 17/6 - 4/(3*A) + 1/(6*pow(A,2) ) );
    return val;
}

double PhotonFlux::getKniehl(double x){
    // similar approximation to Drees, Zeppenfeld, also includes magnetic form factor
    // https://doi.org/10.1016/0370-2693(91)90432-P


    vector<double> c = KniehlConstantsC;
    double y = 1/2 - 2/x +2/pow(x,2);
    double a = 4*pow(M_p,2) /0.71;   // 4.96
    double z = 1 + a/4*(pow(x,2)/(1-x));

    double val = alpha/(2*M_PI) *x* (c[0] *y *log(1 + c[1]/z) -(y+ c[2])*log(1-1/z) + c[3]/(z-1) + (c[4]*y+c[5])/z + (c[6]*y - c[7])/pow(z,2) + (c[8]*y + c[9])/pow(z,3) );
    return val; 

}

double PhotonFlux::getJackson(double x){
    // the analytical solution to the photon flux, not the original source, but the calculation is performed in text book 
    // Classical Electrodynamics b J.D.Jackson 

    double Y = x*M_p *b_min;
    double K0 = gsl_sf_bessel_K0(Y);
    double K1 = gsl_sf_bessel_K1(Y);

    double val = alpha/(M_PI*x) * (2*Y*K0*K1 - pow(Y,2) * ( pow(K1, 2) - pow(K0, 2) ) ); 
    return val;
}

double PhotonFlux::getJackson(double x, double b){
    // same as previous Jackson, although before integration over db^2

    double gamma = getLorentzFactor();
    double omega = omegaFromX(x);
    double K0 = gsl_sf_bessel_K0(b*omega/gamma);
    double K1 = gsl_sf_bessel_K1(b*omega/gamma);
    double val = alpha*pow( omega ,2)/(pow(M_PI*gamma, 2) )*(pow(K1,2) + pow(K0/gamma,2));
    return val;
}

