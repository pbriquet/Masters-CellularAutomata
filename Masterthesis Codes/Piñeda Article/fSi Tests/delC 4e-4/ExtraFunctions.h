//---------------------------------------------------------------------------
#ifndef ExtraFunctionsH
#define ExtraFunctionsH
//---------------------------------------------------------------------------

#include "GlobalDefinitions.h"

double fSScheil(double T);
void getTfar();

double esidtGeral(struct ClusterMatrix *c);

double BirdDel_dl(double vfar, double de, double theta);

double RanzMarshallLenght(double Reynolds, double Prandtl, double Radius);

double delC_Stag(double Re, double Rf, double Reynolds, double Schmidt);

double NewDelDL(double Re, double Rf, double Pe, double delC);

double proj(double v1[2], double v2[2]);

double modulo(double v[2]);

double sumScalarProduct( double coefficient[ 1 + 2*dimmax ], double vector[ 1 + 2*dimmax ], int N);

double IVANTSOV(double Peclet);

double expint( double Peclet );

double Factorial(double n);

double ErrorFunction(double Z);

double LldSphe(double Ri, double Rb, double Pe);

double SoluteDiffusionLength (double El,double velocity,double Ne);

double Tfar();

double Hvar();

double Tboundary( double Tp, double Tfarfield, double Keq, double Hconvection, double dy);

void CountGrainsVF();

void InitiateMatrix();


void RefreshTime();

void DiffusionLenght();

void DelAdjustmentFunction();



#endif
