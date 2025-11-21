#ifndef DipoleCrossSection_h
#define DipoleCrossSection_h


#include "Libraries.h"

using namespace std;


class DipoleCrossSection {
    public:
        DipoleCrossSection();
        ~DipoleCrossSection();


        double xFromS(double s) { return (pow(m_VM, 2) + QSquared )/ s; }
        double sFromX(double x) { return (pow(m_VM, 2) + QSquared )/ x; }
        // functions for different models
        double getGBWOld(double r, double x);
        double getGBWNew(double r, double x);
        double getGBWNewNew(double r, double x);
        double getKST(double r, double x);
        double getGBWDglap( double r, double x); 
        void getData(GridData & s, double* xx, double* yy);
        GridData loadGrid(const string gridName, int xN = 2, int yN = 1);
    private:
        double squareRootOfS = 510; //GeV
        double m_VM = 3.1;  //GeV
        double QSquared = 0;

        double GEV2FM = 5.0677302; // GEV to fm-1 >>   1GeV = 5.0633fm-1  // 1fm = 5.0633GeV-1   ****   1/GEV2FM>   1fm-1=0.197GeV
        double GEV2MB = 0.3893;   // GEV-2 to mb >>   1GeV-2 = 0.3893mb  // mb-1 = 0.3893 GeV2


        // variables for dglap GBW
        GridData mGridData; // data for DGLAP cross section
};

#endif