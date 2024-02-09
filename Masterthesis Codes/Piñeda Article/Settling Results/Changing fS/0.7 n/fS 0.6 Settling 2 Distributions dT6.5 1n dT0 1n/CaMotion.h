//---------------------------------------------------------------------------
#ifndef CaMotionH
#define CaMotionH
//---------------------------------------------------------------------------



double SoluteDiffusionConvecLength( double El, double RelativeVl, double Ne );



double fKvpartition(double el, double beta, double betad);

double fbetad(double esi, double phie, double Senv, double Ssol);

double fCpphie(double phie, double el);

double fbetal(double betad, double el, double phie);

double fbeta(double betad, double betal, double el);
void CalculateAverageVs();

void CalculateAverageVfVdVl();




double EquivalentD( double area );

double DimensionlessDiameter(double dsph, double clusterdensity);

double TerminalVelocity( double DimensionlessV, double clusterdensity);

double DimensionlessVelocity( double DimensionlessD);


double CalculateVelocity( double area , double clusterdensity );

void ClusterFormation();

void ClusterVelocity();


void ClusterMovement();


void CellTranslation();

void VFTranslation();




void NewCdafterConvection();


void printLevenlspiel();

#endif
