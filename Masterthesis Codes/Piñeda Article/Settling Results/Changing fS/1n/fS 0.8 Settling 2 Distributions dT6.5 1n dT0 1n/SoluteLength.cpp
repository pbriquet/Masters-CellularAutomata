//---------------------------------------------------------------------------
#pragma hdrstop

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "SoluteLength.h"
//---------------------------------------------------------------------------
#pragma package(smart_init)

#include "GlobalDefinitions.h"
#include "SetupDatabase.h"
#include "GlobalVariables.h"
#include "ExtraFunctions.h"

double delC(double Re, double Rf, double Rey, double Sch){

        double n;
        double Sh;

        if( Rey == 0.0){
                n = 4.64/14.04;
        }
        else{
            n = (2.0 + 4.64*pow(Rey , -0.24) ) / (3.0 + 14.04 * pow(Rey, -0.24));
        }

        Sh = 2.0*Rf/(Rf - Re) + 2.0/3.0/(1 - pow(Re/Rf, 3.0) )*pow(Rey , n)*pow(Sch , 1.0/3.0);

        return 2.0*Re/Sh;
}

double Clfar(double Cd, double Clave, double Pe, double Re, double Rf, double del){

        double Re_delC = Re + del;
        double Pe2_ = Pe * Re_delC / Re;
        double IvPe = IVANTSOV(Pe);
        double IvPe2 = IVANTSOV(Pe2_);

        double constint = ( pow(Re_delC, 3.0) - pow(Re, 3.0))/3.0*(exp(-Pe) / Pe - IvPe/exp(Pe)/Pe);

        double IntExp = pow( Re/Pe , 3.0)*(Pe2_*exp(-Pe2_) - Pe*exp(-Pe) + exp(-Pe2_) - exp(-Pe));
    
        double IntE1 = pow(Re/Pe , 3.0)*(IvPe2/exp(Pe2_)/Pe2_*pow(Pe2_ , 3.0)/3.0 - IvPe/exp(Pe)/Pe* pow(Pe , 3.0)/3.0 - 1.0/3.0*(pow(Pe2_ , 2.0)*exp(-Pe2_) - pow(Pe, 2.0)*exp(-Pe)) - 2.0/3.0*(Pe2_*exp(-Pe2_) - Pe*exp(-Pe)) - 2.0/3.0*(exp(-Pe2_) - exp(-Pe)));
    
        double A2 = (pow(Re_delC, 3.0) - pow(Re, 3.0)) / (pow(Rf, 3.0) - pow(Re, 3.0));

        double A3 = (pow(Rf, 3.0) - pow(Re_delC, 3.0)) / (pow(Rf, 3.0) - pow(Re, 3.0));
     
        double A1_ = 1.0/(IvPe2/exp(Pe2_)/Re_delC - IvPe/exp(Pe)/Re - exp(-Pe2_)/Re_delC + exp(-Pe)/Re);
     
        double B1 = 3.0*A1_/(pow(Rf , 3.0) - pow(Re , 3.0))*(Pe/Re)*(constint + IntE1 + IntExp);

        return (Clave - A2*Cd + B1*Cd)/(B1 + A3);
}


double Lldave(double Cd, double Clave, double Pe, double Re, double Rf, double Rey, double Sch){

        double CLfar;

        double A1,A2,A3;
        double Pe2, Re_delC;
        double IvPe, IvPe2;

        double Lld;

        Re_delC = Re + delC(Re,Rf,Rey,Sch);

        IvPe = IVANTSOV(Pe);
        IvPe2 = IVANTSOV(Pe*Re_delC/Re);
        Pe2 = Pe*Re_delC/Re;

        Lld = pow(Re , 2.0)*exp(Pe)*(IvPe2/(Re_delC)/exp(Pe2) - IvPe/Re/exp(Pe) - exp(-Pe2)/Re_delC + exp(-Pe)/Re);

        CLfar = Clfar(Cd, Clave, Pe, Re, Rf, Re_delC - Re);

        return Lld*(Cd - Clave)/(Cd - CLfar);




}
