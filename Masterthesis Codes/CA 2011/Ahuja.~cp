#include <math.h>
#include <stdio.h>
#include <stdlib.h>


//---------------------------------------------------------------------------
#pragma hdrstop


#include "Ahuja.h"
#include "GlobalDefinitions.h"
#include "GlobalVariables.h"
//---------------------------------------------------------------------------



#pragma package(smart_init)


double CgAhuja(double phie){

                return 1.2376*log10(phie/0.1556);
}

double CfAhuja(double beta){
                return (2.0*beta*beta*(1.0 - tanh(beta)/beta))/(2.0*beta*beta + 3.0*(1.0 - tanh(beta)/beta));
}

double CdAhuja(double Re, double phie, double beta){

                return 24.0*CfAhuja(beta)/CgAhuja(phie)/Re*(1.0 + CgAhuja(phie)*C1Ahuja*pow(Re, C2Ahuja) ) + ( C3Ahuja/(1.0 + C4Ahuja/Re) );

}

double dCddReAhuja(double Re, double phie, double beta){

                double dRe = 0.00000000000000001;
                double dCD = (CdAhuja(Re + dRe, phie, beta) - CdAhuja(Re,phie,beta))/dRe;
                //double result1 = -24.0*CfAhuja(beta)/CgAhuja(phie)/Re/Re + 24.0*CfAhuja(beta)*C1Ahuja*(C2Ahuja - 1.0)*pow(Re, C2Ahuja - 2.0) + C3Ahuja*C4Ahuja/Re/Re/pow(1.0 + C4Ahuja/Re , 2.0);
                //double result = 24*CfAhuja(beta)/CgAhuja(phie)*(-1.0/Re/Re*(1.0 + CgAhuja(phie)*C1Ahuja*pow(Re,C2Ahuja)))*(CgAhuja(phie)*C1Ahuja*C2Ahuja*pow(Re,C2Ahuja - 1.0))/Re + C3Ahuja*(   1.0/(1.0 + C4Ahuja/(Re + dRe)) - 1.0/(1.0 + C4Ahuja/Re) )/dRe;

                return dCD;
}

double phiReAhuja(double Re, double phie, double beta, double Ve, double de){
                /*
                double dRe = 0.0000001;
                double fRe = Re*Re*CdAhuja(Re,phie,beta) - 8.0*liquiddensity/M_PI*gravity*(soliddensity - liquiddensity)*Ve*clusters[ii].fS;
                double df = ((Re + dRe)*(Re + dRe)*CdAhuja(Re + dRe, phie,beta) - 8.0*liquiddensity/M_PI*gravity*(soliddensity - liquiddensity)*Ve - fRe)/dRe;
                double dfdRe = 2.0*Re*CdAhuja(Re,phie,beta) + Re*Re*dCddReAhuja(Re,phie,beta);
                */
                //double dRe = 0.0000001;
                double fRe = viscosity/liquiddensity*Re/de - sqrt((2.0*clusters[ii].fS*(soliddensity - liquiddensity)*Ve/liquiddensity/(M_PI*de*de/4.0)/CdAhuja(Re,phie,beta)));
                //double df = ((Re + dRe)*(Re + dRe)*CdAhuja(Re + dRe, phie,beta) - 8.0*liquiddensity/M_PI*gravity*(soliddensity - liquiddensity)*Ve - fRe)/dRe;
                double dfdRe = viscosity/liquiddensity/de - sqrt( 2.0*clusters[ii].fS*(soliddensity - liquiddensity)*Ve/liquiddensity/(M_PI*de*de/4.0) )*(1.0/2.0)/pow(CdAhuja(Re,phie,beta), 3.0/2.0);

                return Re - fRe/dfdRe ;
}

double CalculateVClusterAhuja(){

                double fRe = 3.0;
                double fRe1;
                double fRe2;

                double phie;
                double beta;
                double Ve;
                double Re;

                double de;
                double Kperm;
                double Sv;

                double VCr;
                double VC;


                int countador = 0;

                phie = phiex;

                de = sqrt(4.0*clusters[ii].size*dxCA*dyCA/M_PI);
                Ve = M_PI*de*de*de/6.0;

                Sv = 6.0/lambda2;


                if(clusters[ii].fS != 0.0){
                                Kperm = pow(1.0 - clusters[ii].fS, 3.0)/5.0/Sv/Sv/clusters[ii].fS/clusters[ii].fS;
                                beta = de/2.0/sqrt(Kperm);
                }
                else
                                beta = betamax;

                // Initial guess for Reynolds
                if(clusters[ii].V[1] != 0.0){   // if we have a guess
                                Re = liquiddensity*fabs(clusters[ii].V[1])*de/viscosity;
                }
                else{
                                Re = 100000.0;      // if we don't
                }

                while(fRe > 1.0e-6){

                                fRe1 = viscosity/liquiddensity*Re/de;
                                fRe2 = sqrt((2.0*clusters[ii].fS*(soliddensity - liquiddensity)*Ve/liquiddensity/(M_PI*de*de/4.0)/CdAhuja(Re,phie,beta)));
                                fRe = fRe1 - fRe2;


                                Re = phiReAhuja(Re,phie,beta,Ve,de);



                                if(Re < 0.0)
                                                Re = -Re;
                                countador = countador + 1;
                                fRe = viscosity/liquiddensity*Re/de - sqrt((2.0*clusters[ii].fS*(soliddensity - liquiddensity)*Ve/liquiddensity/(M_PI*de*de/4.0)/CdAhuja(Re,phie,beta)));
                                

                }


                VCr = viscosity*Re/liquiddensity/de;

                VC = VCr;

                return VC ;


}

void printAhujaFile(FILE * AAAhujaFile){
                double Re;


                AAAhujaFile = fopen("AAAhujaRes.dat", "w");

                fprintf(AAAhujaFile, "Re;CD;");

                for( Re = 0.1; Re < 1.0e5 ; Re += Re){
                                fprintf(AAAhujaFile, "\n%.10f;%.10f", log10(Re), log10(CdAhuja(Re, 1.0, betamax) ) );
                }

                fclose(AAAhujaFile);

}






