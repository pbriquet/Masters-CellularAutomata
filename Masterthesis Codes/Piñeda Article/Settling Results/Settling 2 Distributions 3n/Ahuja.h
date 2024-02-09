//---------------------------------------------------------------------------
#ifndef AhujaH
#define AhujaH
//---------------------------------------------------------------------------



#define betamax 8.0e8
#define betamin 8.0e-8


#define SphericalAhuja

#define lambda2 6e-6

#ifdef SphericalAhuja
#define phiex 1.0
#endif
#ifdef OctahedronAhuja
#define phiex 0.847
#endif

#define C1Ahuja exp(2.3288-6.4581*phiex + 2.4486*phie*phiex)
#define C2Ahuja (0.0964 + 0.5565*phiex)
#define C3Ahuja exp(4.905 - 13.8944*phiex + 18.4222*phiex*phiex - 10.2599*phiex*phiex*phiex)
#define C4Ahuja exp(1.4681 + 12.2584*phiex - 20.7322*phiex*phiex + 15.8855*phiex*phiex*phiex)




double CgAhuja(double phie);


double CfAhuja(double beta);

double CdAhuja(double Re, double phie, double beta);

double dCddReAhuja(double Re, double phie, double beta);

double phiReAhuja(double Re, double phie, double beta, double Ve, double de);

double CalculateVClusterAhuja();

void printAhujaFile(FILE * AAAhujaFile);



#endif
