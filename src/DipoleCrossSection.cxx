#include "DipoleCrossSection.h"

DipoleCrossSection::DipoleCrossSection(){
  // Constructor implementation if needed
  // Initialize any member variables or resources here

}

DipoleCrossSection::~DipoleCrossSection() {}



double DipoleCrossSection::getGBWOld(double r, double x) {

  // according to arXiv:1901.02664v1 [hep-ph] 

  double sigma0 = 23.03;  //mb
  double x0 = 3.04e-4;
  double lambda = 0.288; //GeV^-2
  double Q0Squared = 1; //GeV^2

  double QSSquared = Q0Squared * pow(x/x0, -lambda);
  return sigma0 * ( 1 - exp(-QSSquared * pow(r, 2) / 4) );
}

double DipoleCrossSection::getGBWNew(double r, double x) {

  // according to arXiv:1901.02664v1 [hep-ph]

  double sigma0 = 23.9;  //mb
  double x0 = 1.11e-4;
  double lambda = 0.287; //GeV^-2
  double Q0Squared = 1; //GeV^2

  double QSSquared = Q0Squared * pow(x/x0, -lambda);
  return sigma0 * ( 1 - exp(-QSSquared * pow(r, 2) / (4 * (1 + QSSquared * pow(r, 2)))) );
}


double DipoleCrossSection::getGBWNewNew(double r, double x){
  // according to code tot_proton.cpp

  double x0= 0.42e-4;
  double lam = 0.248;
  double R02 = 1.* pow(x0/x,lam) * pow(GEV2FM,2) *pow(1.-x,5.6);
  double sigma0 = 27.32 *0.1; // in fm^2


  // standard formula     
  double res = sigma0*(1. - exp(-r*r*R02/4.)); // r in fm, res in fm^2
  //    double res = sigma0*(1. - exp(-r*r/R02)); 
  return res;


}

double DipoleCrossSection::getKST(double r, double x) // r, sq in fm, Gev^2
{
  double s = sFromX(x);
  double s0 = 1000.; // GeV^2
  double R02 = pow(0.88 * pow(s/s0, -0.14),2); // fm^2
  //double B0 = 6.; //GeV-2
  //double aP = 0.25; //GeV-2
  //double mu2 = 1; //GeV2
  double rpi2 = 0.44; // fm^2 
  //double Bel = B0 + 2.*aP*log(s/mu2);
  //double RN2 = Bel - 0.25*R02 - (1./3.)*rpi2; 
  double sigPiPtot = 23.6 * pow(s/s0,0.08) * 0.1 + 1.432 * pow(s/s0, -0.45) * 0.1; // mb -> fm^2
  double sigma0 = sigPiPtot*(1. + 3.*R02/rpi2/8.);
  double res = sigma0*(1. - exp(-r*r/R02));
  // C(s)*r^2 form
  //    double res1 = sigma0*1./R02 * r*r;

  return res;
 }


double DipoleCrossSection::getGBWDglap( double r, double x) // r in fm
{

// K. Golec-Biernat, S. Sapeta, JHEP03 (2018), 102, arXiv:1711.11360

//Limits: x in [1e-6, 1e-2], r in [1e-3, 1e4] in GeV^-1, Q2 in [0.065, 700] in GeV^2.
 //   cout << "input " << r << "\t" << x << endl;

  if (x > 0.01) {
  x = 0.01;
  }

  if (r*GEV2FM < 0.001) {
    r = 0.001/GEV2FM;
  }

  
  double res1[1];
  double xp[2];
  xp[0] = x;
  xp[1] = r*GEV2FM;

  if(mGridData.isData == 1) {
    getData(mGridData, xp, res1);
  }else {
    std::cerr << "ERROR: no data loaded - need to call loadGrid()" << std::endl;
    exit(1);
  }
  
  
  double sigma0 = 22.93; // mb

  double res = res1[0]*sigma0 * 0.1; // mb -> fm^2

//    cout << "output " <<r << "\t" << x << "\t" << res << endl;

    return res;
}

GridData DipoleCrossSection::loadGrid(const string gridName, int xN, int yN) {

  GridData s;
  s.xDim = xN;
  s.yDim = yN;
  s.file = gridName;

  //...open file
  std::ifstream f(s.file.c_str());
  if (!f.is_open()) { // if cannot open
    std::cerr << "ERROR: cannot open file: " << s.file << std::endl;
    exit(1);
  }


  // counting of nonempty lines
  std::string line;
  while (getline(f, line))
    if (!line.empty())
      s.N++;

  //...initialise some arrays arrays
  s.y = new double*[s.yDim]; // 2D
  for (int i = 0; i < s.yDim; ++i) {
    s.y[i] = new double[s.N];
  }

  s.xNum = new int[s.xDim]; // 1D
  s.xLimits = new double*[s.xDim]; // 2D
  for (int i = 0; i < s.xDim; ++i) {
    s.xLimits[i] = new double[2];
  }

  //...temporary arrays
  // will keep all x values
  double **t = new double*[s.xDim]; //2D
  for (int i = 0; i < s.xDim; ++i) {
    t[i] = new double[s.N];
  }
  // will keep y values
  double **ty = new double*[s.yDim]; //2D
  for (int i = 0; i < s.yDim; ++i) {
    ty[i] = new double[s.N];
  }
  // will keep position in array ty
  int **ti = new int*[s.N]; //2D
  for (int i = 0; i < s.N; ++i) {
    ti[i] = new int[s.xDim];
  }



  //...read data
  f.clear(); // retrun to the begining of the file
  f.seekg(0, std::ios::beg); // retrun to the begining of the file
  double* xx = new double[s.xDim];
  for (int i = 0; i < s.N; i++)
  {
    for (int j = 0; j < s.xDim; j++) {
      f >> t[j][i];
    }
    for (int j = 0; j < s.yDim; j++) {
      f >> ty[j][i];
    }
  }

  f.close();

  //...counting number of each x value
  int tmpx;
  int xTot = 0;
  for (int i = 0; i < s.xDim; i++)
  {
    tmpx = 0;
    double xx;
    for (int j = 0; j < s.N; j++) {
      if (j == 0) {
        xx = t[i][j];
        tmpx++;
      }
      else {
        if (t[i][j] > xx) {
          xx = t[i][j];
          tmpx++;
        }
      }
    }
    s.xNum[i] = tmpx;
    xTot += tmpx;
  }

  // this array keep number of xvalues from more cols
  // necessary for array position 
  int* ni = new int[s.xDim];
  for (int i = s.xDim - 1; i >= 0; i--) {
    ni[i] = 1;
    for (int j = s.xDim - 1; j > i; j--) {
      ni[i] *= s.xNum[j];
    }
  }

  // calculation of xND array of ty
  for (int i = 0; i < s.xDim; i++) {
    int n = 0;
    int m = 0;
    for (int j = 0; j < s.N; j++)
    {
      ti[j][i] = m;
      n++;
      if (n > ni[i] - 1) {
        m++;
        n = 0;
      }

      if (m >= s.xNum[i]) {
        m = 0;
        n = 0;
      }
    }
  }


  //...define x values array
  s.x = new double[xTot]; // 1D
  // fill the x values array
  // and find the limits
  int tx = 0;
  for (int i = 0; i < s.xDim; i++)
  {
    double xx;
    for (int j = 0; j < s.N; j++) {
      if (j == 0) {
        xx = t[i][j];
        s.x[tx] = xx;
        tx++;
        s.xLimits[i][0] = xx;
      }
      else {
        if (t[i][j] > xx) {
          xx = t[i][j];
          s.x[tx] = xx;
          tx++;
        }
      }
    }
    s.xLimits[i][1] = xx;
  }

  // next temporary array for changing position
  for (int i = 0; i < s.xDim; i++) {
    ni[i] = 1;
    for (int j = 0; j < i; j++) {
      ni[i] *= s.xNum[j];
    }
  }

  // reposition of values of y for right using with function get
  for (int i = 0; i < s.N; i++) {
    int pos = 0;
    for (int j = 0; j < s.xDim; j++) {
      pos += ti[i][j] * ni[j];
    }
    for (int j = 0; j < s.yDim; j++) {
      s.y[j][pos] = ty[j][i];
    }

  }


  std::cout << "Data loaded:" << std::endl;
  std::cout << "..file: " << s.file << std::endl;
  std::cout << "..number of X values: " << s.xDim << std::endl;
  std::cout << "..number of Y values: " << s.yDim << std::endl;
  std::cout << "..number of lines: " << s.N << std::endl;
  std::cout << "..x values: (i, number, lower limit, upper limit)" << std::endl;


  // data loaded flag
  s.isData = 1;

  //...clean temporary arrs
  delete[] t;
  delete[] ty;
  delete[] ti;
  delete[] ni;

  mGridData = s; // save data for dglap GBW
  return s;

}



void DipoleCrossSection::getData(GridData &s, double* xx, double* yy){
  
  
  
  if (s.isData != 1) {
    std::cerr << "ERROR: no data loaded - need to call ini()" << std::endl;
    return;
  }

  //... check x values and set min/max
  
  for (int i = 0; i < s.xDim; i++) {
    if (xx[i] < s.xLimits[i][0])
      xx[i] = s.xLimits[i][0];
    else if (xx[i] > s.xLimits[i][1])
      xx[i] = s.xLimits[i][1];
  }

  for (int i = 0; i < s.yDim; i++) {
    yy[i] = 0.;
  }



  int* NCOMB = new int[s.xDim];
  double* D = new double[s.xDim];
  int* IENT = new int[s.xDim];
  int KD, M, J, JR, JA, JB;

  KD = 1;
  M = 1;
  JA = 1;

  for (int I = 1; I <= s.xDim; I++)
  {
    NCOMB[I - 1] = 1;
    JB = JA - 1 + s.xNum[I - 1];

    for (J = JA; J <= JB; J++) {
      if (xx[I - 1] <= s.x[J - 1])
        goto g3;
      continue;
    }
    J = JB;

    g3:
      if (J != JA)
        goto g4;

      J = J + 1;

    g4:
      JR = J - 1;
      D[I - 1] = (s.x[J - 1] - xx[I - 1]) / (s.x[J - 1] - s.x[JR - 1]);
      IENT[I - 1] = J - JA;
      KD = KD + IENT[I - 1] * M;
      M = M * s.xNum[I - 1];

    g5:
      JA = JB + 1;

  }

  double FAC;
  int IADR, IFADR, IL;

  g10: 
  FAC = 1.;
  IADR = KD;
  IFADR = 1;
  for (int I = 1; I <= s.xDim; I++)
  {
    if (NCOMB[I - 1] == 0)
      goto g12;

    FAC = FAC*(1. - D[I - 1]);
    goto g15;

    g12:
    FAC = FAC*D[I - 1];
    IADR = IADR - IFADR;

    g15:
    IFADR = IFADR*s.xNum[I - 1];
  }



  for (int k = 0; k < s.yDim; k++) {
    yy[k] = yy[k] + FAC*s.y[k][IADR - 1];
  }
  IL = s.xDim;


  g40: // GO TO 40
  if (NCOMB[IL - 1] == 0)
    goto g80;
  NCOMB[IL - 1] = 0;

  if (IL == s.xDim)
    goto g10;
  IL = IL + 1;

  for (int K = IL; K <= s.xDim; K++) {
    NCOMB[K - 1] = 1;
  }
  goto g10;

  g80: // GO TO 80
  IL = IL - 1;

  if (IL != 0)
    goto g40;
  

  // cleaning...
  delete[] NCOMB;
  delete[] D;
  delete[] IENT;

  //std::cout << std::endl;


}
  