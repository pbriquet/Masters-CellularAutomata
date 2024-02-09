// Deterministic.cpp

#ifndef Deterministic
#define Deterministic

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define AlSi

#define H 250.0
#define Tfar


#ifdef AlSi
#define Co 5.0

#if Co == 5.0
        #define kpart 0.117
        #define pCp 2.35e6 //  J/m³/K
        #define pLf 9.5e8 // J/m³
        #define Lf 372.0e3
        #define Cp 921.5
        #define ml -7.7
        #define Tf 933.0
        #define Tliq Tf + ml*Co
        #define Teut 850.0
        #define Gibbs 0.9e-7
        #define Stability 1.0/4.0/M_PI/M_PI

        #define DiffL 3.0e-9
#endif


#endif


#define Superheat 0.0
#define To Tliq + Superheat

#define Fmax 0.9999

#define Cooling -45.0

#define dtmax 0.001
#define ti 0.0
#define tmax 8.0

// Mesh Properties

#define Lx 0.0079
#define Ly 0.0079
#define VFMeshx 1
#define VFMeshy 1
#define dxVF Lx/VFMeshx
#define dyVF Ly/VFMeshy

//Nucleation Properties

#define dTn 0.0
#define dTsigma 0.0
#define dTmax 0.0
//#define no 2.4e8
#define no nv
#define Rc 1.0e-8


// Abbreviations

#define RUNVFLINE
#define RUNTIME t = ti; t < tmax ; t+=dt

#define CADEResultsTimeInterval 0.01
#define ifCADEResultsTimeInterval (t >= 0.0 && t < 0.0 + dt) || t*(1.0/CADEResultsTimeInterval) - (int) (t*(1.0/CADEResultsTimeInterval) ) >= 0.0 && t*(1.0/CADEResultsTimeInterval) - (int) (t*(1.0/CADEResultsTimeInterval) ) < dt*(1.0/CADEResultsTimeInterval)


//------------------------------------------------------------------



double na = 30.0/Lx/Ly;
double nv = sqrt(M_PI/6)*pow( na, 3.0/2.0);


double Factorial(double n)
{
 int i;
 double w ;

 if (n==0)
   return 1;

 w = 1;
 for( i = 1; i <= n; i++)
   w = w*i;

 return w;
}
//---------------------------------------------------------------------



//---------------------------------------------------------------------
double ErrorFunction(double Z)
{
 int Neg;
 double auxErrorFunc,Lim=1.0E-8,Termo,Sum;
 int n;


 if (Z < 0)
 {
   Neg = 1;
   Z = -Z;
 }
 else
  Neg = 2;


 if (Z > 4)
   auxErrorFunc = 1;
 else
  {
   Sum = Z;
   n = 0;

   do
   {
    n += 1;
    Termo = pow(-1,n)*pow(Z,(2.0*n + 1.0))/(Factorial(n)*(2.0*n + 1.0));
    Sum += Termo;
   }
   while (fabs(Termo) > fabs(Lim * Sum));

   auxErrorFunc = Sum*2.0/pow(M_PI,0.5);
   }

 if (Neg == 1)
   auxErrorFunc = -auxErrorFunc;


 return (auxErrorFunc);
}


double modulo(double Z) {

double aux_modulo;

if (Z<0)
aux_modulo = -Z;
else
aux_modulo = Z;

return (aux_modulo);
}


double IVANTSOV(double Peclet)
{
 double ExpoInt,numerador,denominador,auxIVANTSOV;

 if ((Peclet >= 0) && (Peclet <= 1))
  {
   ExpoInt = -0.57721566 + 0.99999193 * Peclet -
	      0.24991055 * pow(Peclet,2) + 0.05519968 * pow(Peclet,3) -
	      0.00976004 * pow(Peclet,4) + 0.00107857 * pow(Peclet,5) -
	      log(Peclet);

   auxIVANTSOV = Peclet * exp(Peclet) * ExpoInt;
  }
 else
  {
   numerador = pow(Peclet,4) + 8.5733287401 * pow(Peclet,3) +
	      18.059016973 * pow(Peclet,2) + 8.6347608925 * Peclet + 0.2677737343;

   denominador = pow(Peclet,4) + 9.5733223454 * pow(Peclet,3) +
	      25.6329561486 * pow(Peclet,2) + 21.0996530827 * Peclet +
	      3.9584969228;

   auxIVANTSOV = numerador / denominador;
  }

 return (auxIVANTSOV);

}



double expint( double Peclet ){
        return -0.57721566  +
		0.99999193*Peclet -
		0.24991055*pow(Peclet,2) +
		0.05519968*pow(Peclet,3) -
		0.00976004*pow(Peclet,4) +
		0.00107857*pow(Peclet,5) -
		log(Peclet); /*ln*/
}




//------------------------------------------------------------------

double LldSphe(double Ri, double Rb, double Pe)
{
 double aux1,aux2,aux3,auxLldSphe;

 aux1 = aux1;
 if (Pe > 1.0E-5)
 {
  aux1 = (((Rb * Ri) / Pe) + pow(Ri/Pe,2) - pow(Rb,2)) * exp(-Pe * ((Rb / Ri) - 1.0));

  aux2 = ((Ri * Ri) / Pe) + pow(Ri/Pe,2) - (pow(Rb,3) / Ri);

  aux3 = (Pe * pow(Rb,3) / Ri) * (exp(-Pe*((Rb / Ri) - 1.0)) *
	  IVANTSOV(Pe * Rb / Ri) / (Pe * Rb / Ri) - IVANTSOV(Pe) / Pe);

  auxLldSphe = ( pow(Ri,2)/(pow(Rb,3)-pow(Ri,3)) ) * (aux1 - aux2 + aux3);
 }

 else if (Pe <= 1.0E-5)
  auxLldSphe = Ri*(1 -(3.0/2.0)*(Ri*(pow(Rb,2)-pow(Ri,2))/(pow(Rb,3)-pow(Ri,3))));

 return (auxLldSphe);

}




//---------------------------------------------------------------------





//-Cálculo de del_e---------------------------------------------------------

double SoluteDiffusionLength (double El,double velocity,double Ne)
{
 double auxPeclet,Rcell,Re,auxDiffusionLentgth,Eg;

 if(velocity < 0.0)
        velocity = velocity;

 Eg = 1 - El;

 if (Eg > 0.97 )
   Eg = 0.97;
 else if (Eg < 1E-12)
   Eg = 1E-12;

 Rcell = pow(3.0/(4.0*M_PI*Ne),1.0/3.0);
 Re = Rcell*pow(Eg,(1.0/3.0));
 auxPeclet = velocity*Re/DiffL;

 auxDiffusionLentgth = LldSphe(Re,Rcell,auxPeclet);

 return (auxDiffusionLentgth);

}

//------------------------------------------------------------------





// Global Variables

double t;
double dt;

double Ve;


struct Matrix{

        double T;
        double Cl;
        double Cd;
        double Cs;
        double Fs;
        double Fd;
        double Fl;
        double Se;
        double ne;


};




// Auxiliar Functions

double InverseIvantsov(double Cd, double Cl){

        double omega;

        omega = (Cd - Cl)/Cd/(1.0 - kpart);

        if(omega < 0.0){
                return 0.0;
        }
        else
                return 0.4567*pow(omega/(1.0 - omega) , 1.195);
}


double Vcol(double Cd, double Cl){

        return (4.0*Stability*DiffL*ml*(kpart - 1.0)*Cd)/Gibbs*pow(InverseIvantsov(Cd,Cl), 2.0);

}











int main()
{

        struct Matrix VF;
        struct Matrix VFn;

        // Initiate Matrix VF and VFn
        VFn.T = VF.T = To;
        VFn.Cl = VF.Cl = Co;
        VFn.Cd = VF.Cd = Co;
        VFn.Cs = VF.Cs = Co;
        VFn.Fl = VF.Fl = 1.0;
        VFn.Fd = VF.Fd = 0.0;
        VFn.Fs = VF.Fs = 0.0;
        VFn.Se = VF.Se = 0.0;
        VFn.ne = VF.ne = 0.0;

        double tdimless;
        double Tdimless;

        double JDld = 0.0;
        double del_ld;
        double Rf = 0.0;

        double errorp = 0.0;

        bool instantnuc = 1;

        FILE *TdlXtdl;
        TdlXtdl = fopen("Txtdim.dat", "w");

        FILE *TFile;
        TFile = fopen("T.dat", "w");

        fprintf(TdlXtdl, "tdim;Tdim;");
        fprintf(TFile, "t;T;") ;


        FILE *CLFile;
        CLFile = fopen("Cl.dat", "w");

        FILE *CdFile;
        CdFile = fopen("Cd.dat", "w");

        FILE *egesFile;
        egesFile = fopen("eges.dat", "w");

        FILE *SeFile;
        SeFile = fopen("Se.dat", "w");


        fprintf(CLFile, "t;Cl;");
        fprintf(CdFile, "t;Cd;");
        fprintf(egesFile, "t;eg;es;");
        fprintf(SeFile, "t;Se;");


        dt = dtmax;

        for(RUNTIME){

                //PrintResults

                if(ifCADEResultsTimeInterval){
                tdimless = -pCp*Cooling*t/pLf;
                Tdimless = (To - VF.T)/(To - Teut);
                fprintf(TdlXtdl, "\n%f;%f;", tdimless, Tdimless);
                fprintf(TFile, "\n%f;%f;", t, VF.T);
                fprintf(CLFile, "\n%f;%f;", t, VF.Cl);
                fprintf(CdFile, "\n%f;%f;", t, VF.Cd);
                fprintf(egesFile, "\n%f;%f;%f", t, 1.0 - VF.Fl, VF.Fs);
                fprintf(SeFile, "\n%f;%f;", t, VF.Se);

                }



                //Calculate Variations in Time
                if(VF.Cl < Co)
                        JDld = JDld;

                Ve = Vcol(VF.Cd, VF.Cl);

                if(VF.ne != 0.0){
                        del_ld = SoluteDiffusionLength (VF.Fl,Ve,VF.ne);
                        if( VF.Cd >= VF.Cl)
                        JDld = VF.Se*DiffL/del_ld*(VF.Cd - VF.Cl)*dt;
                        else
                        JDld = 0.0;
                } // End of if(VF.Ne != 0.0)
                else{ // if(VF.Ne == 0.0)
                        JDld = 0.0;
                } // End of else if(VF.Ne == 0.0)

                //if(VF.T > Teut){

                        VFn.Fs = VF.Fs + ( VF.Fd*(VF.T + Cooling*dt - Tf - ml*VF.Cd) + ml*JDld)/(ml*(1.0 - kpart)*VF.Cd - VF.Fd*Lf/Cp );

                        if( VFn.Fs < 0.0){
                                VFn.Fs = 0.0;
                        }
                        if( VFn.Fs > Fmax){
                                VFn.Fs = Fmax;
                        }

                        VFn.T = VF.T + Cooling*dt + Lf/Cp*(VFn.Fs - VF.Fs);

                        errorp = -Cooling*dt*Cp/Lf - (VFn.Fs - VF.Fs);

                        if(VFn.T > VF.T)
                                VF.T = VF.T;

                        VFn.Fl = VF.Fl - VF.Se*Ve*dt;

                        if(VFn.Fl < 1.0 - Fmax){
                                VFn.Fl = 1.0 - Fmax;
                        }


                        #if dTsigma != 0.0

                        #endif

                        #if dTsigma == 0.0
                                if(instantnuc == 1 && (Tliq - VFn.T) > 0.0 && (Tliq - VFn.T) >= dTn){
                                        VFn.ne = VF.ne + instantnuc*no;
                                        instantnuc = 0;
                                }
                        #endif


                        if(VF.ne == 0 && VFn.ne != 0.0) {
                                Rf = pow(3.0/4.0/M_PI/VFn.ne , 1.0/3.0);
                                VFn.Se = VFn.ne*Rc*Rc*4.0*M_PI;
                                //VFn.Fl = 1.0 - pow(Rf*VFn.Se/3.0 , 3.0/2.0);

                        }
                        if(VFn.ne != 0.0 && VF.ne != 0.0){
                                Rf = pow(3.0/4.0/M_PI/VFn.ne , 1.0/3.0);
                                VFn.Se = 3.0*pow(1.0 - VFn.Fl, 2.0/3.0)/Rf ;
                        }



                        VFn.Fd = 1.0 - VFn.Fl - VFn.Fs;

                        if(VFn.Fd < 0.0){
                                VFn.Fd = 0.0;
                                VFn.Fl = 1.0 - VFn.Fs - VFn.Fd;
                        }
                        if(VFn.Fd > Fmax){
                                VFn.Fd = Fmax;
                                VFn.Fl = 1.0 - VFn.Fs - VFn.Fd;
                        }






                        if(VF.Fd != 0.0)
                                //VFn.Cd = VF.Cd + ( (1.0 - kpart)*VF.Cd*(VFn.Fs - VF.Fs) - JDld )/VF.Fd;
                                VFn.Cd = (VFn.T - Tf)/ml;
                        else
                                VFn.Cd = (VFn.T - Tf)/ml;

                        if(VFn.Fl != 0.0)
                                VFn.Cl = (VF.Fl*VF.Cl + VF.Cd*(VFn.Fl - VF.Fl) + JDld)/VFn.Fl;








                //} // End of if(VF.T > Teut)
                else{ // VF.T <=Teut

                } // End of else if(VF.T <= Teut)







                //RefreshTime

                VF.T = VFn.T;
                VF.Cl = VFn.Cl;
                VF.Cd = VFn.Cd;
                VF.Cs = VFn.Cs;
                VF.Fl = VFn.Fl;
                VF.Fd = VFn.Fd;
                VF.Fs = VFn.Fs;
                VF.Se = VFn.Se;
                VF.ne = VFn.ne;

        }

        return 0;
}


#endif
