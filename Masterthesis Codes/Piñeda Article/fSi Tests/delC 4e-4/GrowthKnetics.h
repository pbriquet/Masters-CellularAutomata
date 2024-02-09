//---------------------------------------------------------------------------
#ifndef GrowthKneticsH
#define GrowthKneticsH
//---------------------------------------------------------------------------


#ifdef eglocal
void CalculateAmax();
#endif

double ThermalVelocity(double dT);

double SolutalVelocity( double Cd, double Cl);

double DendriticVelocity();


void CAGrowth();





#endif
