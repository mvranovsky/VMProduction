#ifndef WaveFunction_h
#define WaveFunction_h


#include "Libraries.h"


using namespace std;



class WaveFunction {
    public:
        WaveFunction(string pathToGrid);
        ~WaveFunction();

        void getData(GridData & s, double* xx, double* yy);
        GridData loadGrid(const string gridName, int xN = 2, int yN = 1);

    private:
        string mGridName;
        GridData mGridData;

};

#endif