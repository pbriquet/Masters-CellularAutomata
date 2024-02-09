
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
//---------------------------------------------------------------------------
#pragma hdrstop

#include "DiffusionLenghts.h"
//---------------------------------------------------------------------------

#pragma package(smart_init)


#include "GlobalDefinitions.h"
#include "GlobalVariables.h"

#include "GrowthKnetics.h"
#include "ExtraFunctions.h"


double BirdDiffusionLength(double El, double vfar, double Ne){
                double Reynolds;
                double Schmidt;
                double de, Eg, Rcell;
                double ftheta, theta;

                Eg = 1 - El;

                if (Eg > 0.97 )
                                Eg = 0.97;
                else if (Eg < 1E-12)
                                Eg = 1E-12;

                Rcell = pow(3.0/(4.0*M_PI*Ne),1.0/3.0);
                de = 2.0*Rcell*pow(Eg,(1.0/3.0));

                Reynolds = liquiddensity*de*vfar/viscosity;
                Schmidt = viscosity/liquiddensity/DiffLiq;

                if(side == WEST){
                                theta = M_PI/2.0;
                                ftheta = pow(M_PI - theta + sin(2.0*theta)/2.0 , 1.0/3.0)/sin(theta);

                }
                if(side == EAST){
                                theta = M_PI/2.0;
                                ftheta = pow(M_PI - theta + sin(2.0*theta)/2.0 , 1.0/3.0)/sin(theta);
                }
                if(side == NORTH){
                                ftheta = 1464.6;
                }
                if(side == SOUTH){
                                ftheta = 0.87343322;
                }

                del_dl = pow(3.0/4.0 , 1.0/3.0)*de*pow( Reynolds*Schmidt , -1.0/3.0)*ftheta;

                return del_dl;

}

void ChosenDel_dl(){

                velocity = DendriticVelocity();

                #ifdef DiffusiveLength
                del_dl = SoluteDiffusionLength( VfP.fL, velocity, VfP.Ne );
                #endif
                #ifdef BirdLength
                del_dl = BirdDiffusionLength(VfP.fL, fabs(clusters[ClusterP].V[1]), VfP.Ne);
                #endif

                #ifdef CriteriaLength
                del_dl = CriteriaDiffusionLength();
                #endif
                #ifdef WangBeckermann1996Legth
                del_dl = SoluteDiffusionConvecLength( VfP.fL, fabs(VfP.VL - VfP.Vs), VfP.Ne );
                #endif

}


