#ifndef Libraries_h
#define Libraries_h


#include <vector>
#include <iostream>
#include <cmath>
#include <string>
#include <fstream>
#include <sstream> 
#include <algorithm> 
#include <stdio.h> 
#include <stdlib.h> 
#include <vector> 
#include <fstream> 
#include <cmath>
#include <math.h> 
#include <cstdlib>
#include <bitset>
#include <map>
#include <memory>
#include <unordered_map>
#include <algorithm>
#include <set>

struct GridData {
    double ** y;        // y values
    double * x;         // x values
    double ** xLimits;  // limits of each x values
    int * xNum;         // number of values of x values
    int isData;         // flag is data loaded
    int xDim;           // number of x values
    int yDim;           // number of y values
    int N;              // number of lines
    std::string file;   // name of file
};

#endif