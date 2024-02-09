//---------------------------------------------------------------------------
#pragma hdrstop

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "ExtraFunctions.h"
//---------------------------------------------------------------------------
#pragma package(smart_init)


#include "GlobalDefinitions.cpp"
#include "SetupDatabase.cpp"
#include "GlobalVariables.h"

#include "OutputFiles.h"
#include "CaMotion.h"

#include "InterpolationFunctions.h"

#include "GrowthKnetics.h"


double fSScheil(double T){
        double fS;
        fS = 1.0 - pow((T-Tf)/(Tliq - Tf), 1.0/(kpart - 1.0));
        return fS;
}

void esidtGeral(struct ClusterMatrix *c){
        c->esidt = c->fS - c->fS*3.0/c->de*c->ve*dt + liquiddensity/soliddensity*(1 - c->fS)/(1 - kpart)/c->Cdt*c->dCddt*dt + liquiddensity/soliddensity*c->Ae/c->Ve*Dl/del_C*(c->Cdt - Co)/(1 - kpart)/c->Cdt*dt;
}

void getTfar(){
        FILE *TfarInput;
        char ah[20], *s;
        int u, v;
        TfarInput = fopen("TfarInput.txt","r");
        for(u = 0; u < 28; u ++){
                for(v = 0; ah[v-1] != '\t'; v++)
                        fscanf(TfarInput, "%c", &ah[v]);

                T4vst[u][0] = strtod(ah, &s);

                for(v = 0; ah[v-1] != '\n'; v++)
                        fscanf(TfarInput, "%c", &ah[v]);

                T4vst[u][1] = strtod(ah, &s);
        }

}

double BirdDel_dl(double vfar, double de, double theta){

           double Schmidt = viscosity/DiffLiq/liquiddensity;
           double Reynolds = liquiddensity*vfar*de/viscosity;

           double result;

           Reynolds = liquiddensity*vfar*de/viscosity;
           Schmidt = viscosity/DiffLiq/liquiddensity;

           result = pow(3.0/4.0 , 1.0/3.0)*de*pow(Reynolds*Schmidt , -1.0/3.0)*0.87343322;

           return result*Reynolds*Schmidt/de;

}

double RanzMarshallLenght(double Reynolds, double Prandtl, double Radius){

        return 2.0*Radius/(2.0 + 0.6*pow(Reynolds, 1.0/2.0)*pow(Prandtl,  1.0/3.0));

}

double delC_Stag(double Re, double Rf, double Reynolds, double Schmidt){

        double nexpo;


        if( Reynolds < 1.0e-8){
                nexpo = 4.64/14.04;
        }
        else{
                nexpo = (2.0 + 4.64 * pow(Reynolds , -0.24) ) / (3.0 + 14.04 * pow(Reynolds , -0.24));
        }

        double Sherwood = 2.0*Rf/(Rf - Re) + 2.0/3.0/(1.0 - pow(Re/Rf , 3.0) ) * pow(Reynolds , nexpo)*pow( Schmidt , 1.0/3.0);

        if(Sherwood == 2.0){
                return (Rf - Re);
        }
        return (2.0*Re/Sherwood)/(1.0 - 2.0/Sherwood);

}


double NewDelDL(double Re, double Rf, double Pe, double delC){

        double RedelC = Re + delC;

        Re = Re;

        if( Pe == 0.0){
                return Re*Re/(Rf*Rf*Rf - Re*Re*Re)*(Rf*Rf*Rf/Re - Rf*Rf*Rf/RedelC - RedelC*RedelC/2.0 + Re*Re/2.0);
        }
        else{
                return Re*Re/(Rf*Rf*Rf - Re*Re*Re)*((Rf*Rf*Rf - RedelC*RedelC*RedelC) * (exp(-Pe*(RedelC/Re - 1.0)) * (IVANTSOV(Pe*RedelC/Re) - 1.0)/RedelC + (1.0 - IVANTSOV(Pe))/Re) + exp(-Pe *(RedelC/Re - 1.0))*(RedelC*Re/Pe + Re*Re/Pe/Pe - RedelC*RedelC) - (Re*Re/Pe + Re*Re/Pe/Pe - RedelC*RedelC*RedelC/Re) + Pe*RedelC*RedelC*RedelC/Re*(IVANTSOV(Pe*RedelC/Re)*exp(-Pe*(RedelC/Re - 1.0))/(Pe*RedelC/Re) - IVANTSOV(Pe)/Pe));
        }



}






double proj(double v1[2], double v2[2]){        /*Retorna a projeção de um vetor v1, sobre um vetor unitário v2 */
        return v1[0]*v2[0] + v1[1]*v2[1] ;
}

double modulo(double v[2]){
        return sqrt( pow( v[0], 2.0) + pow( v[1], 2.0) );
}

double sumScalarProduct( double coefficient[ 1 + 2*dimmax ], double vector[ 1 + 2*dimmax ], int N){

        double sum = 0;
        for( int k = 0; k < N; k++){
                sum += coefficient[k]*vector[k];
        }

        return sum;

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



//------------------------------------------------------------------

double LldSphe(double Ri, double Rb, double Pe)
{
 double aux1,aux2,aux3,auxLldSphe;

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

Ne = Ne/mmiy/mmix;
#ifdef SeCACorrectionDiagonalSolutalDiffusionLenght
velocity = velocity/sqrt(2.0);
#endif
 //if(velocity != 0.0){
 double auxPeclet,Rcell,Re,auxDiffusionLentgth,Eg,Velocity,result;
 Velocity = velocity;
                if(Ne == 0.0)
                                return 1.0;

 /*
 if(ii != my - 1 && CaS.state == 0)
        Velocity = Velocity + ( fabs(clusters[ClusterP].V[1]) + fabs(VfP.Vd) );
 if(ii != 0 && CaN.state == 0)
        Velocity = Velocity + ( fabs(clusters[ClusterP].V[1]) + fabs(VfP.Vd) );

 if(Velocity < 0.0)
        ii=ii;
 */

 Eg = 1 - El;

 if (Eg > 0.97 )
   Eg = 0.97;
 else if (Eg < 1E-12)
   Eg = 1E-12;

 Rcell = pow(3.0/(4.0*M_PI*Ne),1.0/3.0);
 Re = Rcell*pow(Eg,(1.0/3.0));
 auxPeclet = Velocity*Re/Dl;

 auxDiffusionLentgth = LldSphe(Re,Rcell,auxPeclet);

 //result = auxDiffusionLentgth*exp(-1500.0*( fabs(clusters[ClusterP].V[1]) + fabs(VfP.Vd) ) );
 result = auxDiffusionLentgth;

 if(result < 0)
        result = result;

 return result;

 //return (1.0e-6*exp(-2000.0*(0.005 + 0.001)));
 //return 1.0e-6;
}



//------------------------------------------------------------------

#define tpause 104.0
#define tpol 109.585


double Tfar(){
        double tconvert = t/60.0 + tpol;
        double Tconvert;
        /*
        if( tconvert < tpol ){
                Tconvert = 679.2036842;
        }
        else{
                Tconvert = 183571.516402954-4421.68668232945*tconvert+35.6086147296172*tconvert*tconvert-0.09571550158173*tconvert*tconvert*tconvert;
        }
        */

        static int upos;

        if(t >= T4vst[0][0] && t < T4vst[1][0])
                upos = 0;

        if(t >= T4vst[upos+1][0])
                upos ++;

        Tconvert = T4vst[upos][1] + (T4vst[upos+1][1] - T4vst[upos][1])/(T4vst[upos+1][0] - T4vst[upos][0])*(t - T4vst[upos][0]);
        return Tconvert + 273.0;
}

double Hvar(){
        double H0 = 600;
        double Hfinal = 1850;
        double deltt = 400.0;
        double Hcalc;
        if(t < deltt)
                Hcalc = H0 + (Hfinal - H0)/deltt*t;
        else
                Hcalc = Hfinal;
        return Hcalc;
}


double Tboundary( double Tp, double Tfarfield, double Keq, double Hconvection, double dy){

        return ( Hconvection*Tfarfield + 2.0*Keq/dy*Tp )/(Hconvection + 2.0*Keq/dy );

        /*
        #ifdef VirtualZeroDimensions
                if( Hconvection != 0.0 )
                                return Tfarfield - CoolingRate*dyVF/Hconvection + Keq/soliddensity/Cp/Hconvection/dyVF*(VF[mmay - 2].T - VF[mmay -1].T);
                else
                                return ( Hconvection*Tfarfield + 2.0*Keq/dy*Tp )/(Hconvection + 2.0*Keq/dy );
        #endif
        */


}

void CountGrainsVF(){

        if( fabs(CaP.state) != 0 && fabs(CaP.state) != stateholder){
                                                                        // If the state of the last grain found is different from the the actual,
                  // NeSearch already keep his position. And then it is verified if he already belonged to the array
                stateholder = fabs(CaP.state);

                NeKeeper[VfP.Ne] = stateholder;

                VfP.Ne++;

                for( int k = 0; k < VfP.Ne - 1; k++){
                        if(NeKeeper[k] == stateholder){
                        VfP.Ne --;
                        NeKeeper[VfP.Ne] = 0;
                        break;
                        }
                }

        }
        if( CaP.dT >= 0.0)
                VfP.Nuc++;


}

















void InitiateMatrix(){

        if((auxpointmacro=malloc(mmay*mmax*sizeof(struct MacroMatrix)))!=NULL){
                VF = (struct MacroMatrix*) auxpointmacro;
        }

        if((auxpointmacroaux=malloc(mmay*mmax*sizeof(struct MacroMatrix)))!=NULL){
                VFaux = (struct MacroMatrix*) auxpointmacroaux;
        }
                // 3.3.4) Registro dos envelopes

        if((auxpointclusters=malloc(mx*my*sizeof(struct ClusterMatrix)))!=NULL){
                clusters = (struct ClusterMatrix*) auxpointclusters;
        }

                 // 3.3.1) Matriz microscópica

        if((auxpointmicro=malloc(my*mx*sizeof(struct MicroMatrix)))!=NULL){
                CA = (struct MicroMatrix*) auxpointmicro;
        }

        if((auxpointmicroaux=malloc(my*mx*sizeof(struct MicroMatrix)))!=NULL){
                CAaux = (struct MicroMatrix*) auxpointmicroaux;
        }



        double Teste = To;

        NTotalGrains = 0;

        for(RunVFLine){      /* Inicializar matriz macroscópica */
                for(RunVFCol){

                        VfP.T = To;
                        VfP.fS = 0.0;
                        VfP.CL = Co;
                        #ifdef GradCInitial
                        VfP.CL = Co + ia;
                        #endif
                        VfP.dH = 0.0;
                        VfP.dT = 0.0;
                        VfP.dfS = 0.0;
                        VfP.CLdt = Co;
                        VfP.VL = 0.0;
                        VfP.state = 0;
                        VfP.Vd = 0.0;
                        VfP.Vs = 0.0;
                        VfP.Vf = 0.0;
                        VfP.Se = 0.0;


                        VfP.Cd = Co;
                        VfP.fL = 1.0;
                        VfP.dfL = 0.0;

                        VfP.Ne = 0.0;
                        VfP.Nuc = 0.0;

                        #ifdef Eutetic
                        VfP.fE = 0.0;
                        VfP.dfE = 0.0;
                        #endif


                }
        }

        for(RunVFLine){
        for(RunVFCol){
        for(RunCAinVFLine){  /* Zerar duas matrizes microscópicas */
                for(RunCAinVFCol){

                        CaP.state = 0;
                        CaP.fS = 0.0;
                        CaP.dfS = 0.0;
                        CaP.theta = 0.0;
                        CaP.L[0] = 0.0;
                        CaP.L[1] = 0.0;
                        CaP.L[2] = 0.0;
                        CaP.L[3] = 0.0;
                        CaP.dC[0] = 0.0;
                        CaP.dC[1] = 0.0;
                        CaP.dT = -1.0;
                        CaP.Cd = VfP.CL;
                        CaP.Cs = kpart*CaP.Cd;
                        CaP.Csdt = kpart*CaP.Cd;

                        #ifdef KeepClinCA
                        CaP.Cl = 0.0;
                        #endif

                        #ifdef eglocal
                        CaP.Amin = 0.0;
                        CaP.Amax = 0.0;
                        #endif

                        #ifdef Eutetic
                        CaP.fE = 0.0;
                        CaP.dfE = 0.0;
                        #endif

                          /* Zerar clusters de envelopes */

                        clusters[ii*mx + ji].thetaclass = 0;
                        clusters[ii*mx + ji].V[0] = 0.0;
                        clusters[ii*mx + ji].V[1] = 0.0;
                        clusters[ii*mx + ji].d[0] = 0.0;
                        clusters[ii*mx + ji].d[1] = 0.0;
                        clusters[ii*mx + ji].left = 0;
                        clusters[ii*mx + ji].right = 0;
                        clusters[ii*mx + ji].up = 0;
                        clusters[ii*mx + ji].down = 0;
                        clusters[ii*mx + ji].fixed = 0;
                        clusters[ii*mx + ji].cluster = ii*mx + ji;
                        clusters[ii*mx + ji].size = 0.0;
                        clusters[ii*mx + ji].fS = 0.0;
                        clusters[ii*mx + ji].Dave = 0.0;

                        #ifdef esiCluster
                        clusters[ii*mx + ji].dCddt = Co;
                        clusters[ii*mx + ji].Cdt = Co;
                        clusters[ii*mx + ji].de = 0.0;
                        clusters[ii*mx + ji].ve = 0.0;
                        clusters[ii*mx + ji].Ae = 0.0;
                        clusters[ii*mx + ji].Ve = 0.0;
                        #endif



                        CaPaux.state = 0;
                        CaPaux.fS = 0.0;
                        CaPaux.dfS = 0.0;
                        CaPaux.theta = 0.0;
                        CaPaux.L[0] = 0.0;
                        CaPaux.L[1] = 0.0;
                        CaPaux.L[2] = 0.0;
                        CaPaux.L[3] = 0.0;
                        CaPaux.dC[0] = 0.0;
                        CaPaux.dC[1] = 0.0;
                        CaPaux.dT = -1.0;
                        CaPaux.Cd = VfP.CL;
                        CaPaux.Cs = kpart*CaP.Cd;
                        CaPaux.Csdt = kpart*CaP.Cd;

                        #ifdef eglocal
                        CaPaux.Amin = 0.0;
                        CaPaux.Amax = 0.0;
                        #endif


                        CaPaux.fE = 0.0;
                        CaPaux.dfE = 0.0;


                }
        }
        }
        }






        for( int k = 0; k < 1 + 2*dimmax; k++){
                aT[k] = 0;
                aflCl[k] = 0;
                afl[k] = 0;
                afdCd[k] = 0;
                afd[k] = 0;

        }

        for( int k = 0; k < mmicrox*mmicroy; k++){
                NeKeeper[k] = 0;
        }



        dt = dtmax;


        #ifdef ALL1GRAIN
                for(RunCALine){
                for(RunCACol){
                        CaP.state = 2;
                        CaP.theta = 0.0;
                        CaP.fS = fSCAConst;
                        CaP.Cd = CdCAConst;
                        clusters[0].size = mx*my;
                        NTotalGrains = 1;
                }
                }
        #endif


        #ifdef INITIALGRAIN

                int countsize = 0;

                initialii = (my/8.0 - sizeofInitialGrain/2) + (sizeofInitialGrain/2);
                //initialii = my/2;
                initialposition = initialii*mx + initialji;
                CA[initialposition].state = 2;
                CA[initialposition].fS = fSCAConst;
                CA[initialposition].Cd = CdCAConst;
                countsize = 1;

                /*
                for(ii = initialii - (int) (sizeofInitialGrain/2) ; ii < initialii + (int) (sizeofInitialGrain/2) ; ii++){
                                for(ji = initialji - (int) (sizeofInitialGrain/2) ; ji < initialji + (int) (sizeofInitialGrain/2) ; ji++){
                                                CaP = CA[initialposition];
                                                countsize ++;
                                }
                }
                */

                for(RunCALine){
                        for(RunCACol){
                                if(pow( ii - initialii ,2.0)*dyCA*dyCA + pow( ji - initialji, 2.0)*dxCA*dxCA < pow(sizeofInitialGrain/2.0 , 2.0)){
                                        CaP = CA[initialposition];
                                        countsize ++;
                                }
                        }
                }


                NTotalGrains = 1;
                clusters[0].size = countsize;
                clusters[0].fS = fSCAConst;
                clusters[0].V[1] = -CalculateVelocity( clusters[0].size*dxCA*dyCA , clusters[0].fS*soliddensity + (1.0 - clusters[0].fS)*liquiddensity );

                VFTranslation();
                for(RunVFLine){
                        for(RunVFCol){
                                stateholder = 0;
                                for(RunCAinVFLine){
                                        for(RunCAinVFCol){
                                                CountGrainsVF();
                                        }
                                }
                        }
                }
                for(RunVFLine){
                                for(RunVFCol){
                                                VfP.fL = VfP.fLCA;
                                }
                }

                #ifdef FlowPartitionCoefficient
                for(RunVFLine){
                        for(RunVFCol){
                                VfP.Vs = 0.0;
                                for(RunCAinVFLine){
                                        for(RunCAinVFCol){

                                                CalculateAverageVs();

                                        }
                                }
                                CalculateAverageVfVdVl();
                        }
                }
                #endif
        #endif


        #ifdef INITIALGRAINS

                int kgrao;
                
                int countsize[NumberofInitialGrainsx*NumberofInitialGrainsy];
                double xo[NumberofInitialGrainsx*NumberofInitialGrainsy], yo[NumberofInitialGrainsx*NumberofInitialGrainsy];
                double dxCell, dyCell, xCA, yCA;
                dxCell = dxVF/NumberofInitialGrainsx;
                //dyCell = dyVF/NumberofInitialGrainsy;
                dyCell = dxVF/NumberofInitialGrainsy;
                for(kgrao=0;kgrao<NumberofInitialGrainsx*NumberofInitialGrainsy;kgrao++){
                                countsize[kgrao]=0;
                                xo[kgrao] = (kgrao%NumberofInitialGrainsx)*dxCell + 0.5*dxCell;
                                yo[kgrao] = Ly - ((kgrao/NumberofInitialGrainsy)*dyCell + 0.5*dyCell);

                }


                for(RunCALine){
                        for(RunCACol){
                                xCA = (ii + 0.5)*dxCA;
                                yCA = (my - 1 - ji + 0.5)*dyCA;
                                for( kgrao=0;kgrao<NumberofInitialGrainsx*NumberofInitialGrainsy;kgrao++){

                                                if(pow( xCA - xo[kgrao] ,2.0) + pow( yCA - yo[kgrao], 2.0) < pow(sizeofInitialGrain/2.0 , 2.0) ){
                                                                CaP.state = kgrao + 2;
                                                                CaP.fS = fSCAConst;
                                                                CaP.Cd = CdCAConst;
                                                                countsize[kgrao]++;
                                                }

                                }
                        }
                }



                NTotalGrains = NumberofInitialGrainsx*NumberofInitialGrainsy;

                for(kgrao=0;kgrao<NumberofInitialGrainsx*NumberofInitialGrainsy;kgrao++){
                                //clusters[kgrao].size = countsize[kgrao];
                                clusters[kgrao].size = (M_PI*pow(Lx*egproj/NumberofInitialGrainsx, 2.0))/(dxCA*dyCA);
                                clusters[kgrao].fS = fSCAConst;
                                clusters[kgrao].V[1] = -CalculateVelocity( clusters[kgrao].size*dxCA*dyCA , clusters[kgrao].fS*soliddensity + (1.0 - clusters[kgrao].fS)*liquiddensity );


                }

                VFTranslation();
                for(RunVFLine){
                        for(RunVFCol){
                                stateholder = 0;
                                for(RunCAinVFLine){
                                        for(RunCAinVFCol){
                                                CountGrainsVF();
                                        }
                                }
                        }
                }
                for(RunVFLine){
                                for(RunVFCol){
                                                VfP.fL = VfP.fLCA;
                                }
                }

        #endif

        // Initiate Files

        InitiateFiles();

        #ifdef MartoranoLiquidBoundaryCondition
        for(RunCACol){
                BoundaryConditionHolder[ji].on = 0;
                BoundaryConditionHolder[ji].edVd = 0.0;                
        }
        #endif



}


    void RefreshTime(){

        for(RunVFLine){
                for(RunVFCol){

                        VfP.T += VfP.dT;
                        VfP.fS += VfP.dfS;
                        VfP.fL += VfP.dfL;
                        VfP.CL = VfP.CLdt;

                                VfP.fE += VfP.dfE;


                        VfP.dT = 0.0;
                        VfP.dfS = 0.0;
                        VfP.dfL = 0.0;
                        VfP.CLdt = VfP.CL;

                        for(RunCAinVFLine){
                                for(RunCAinVFCol){

                                        CaP.fS += CaP.dfS;
                                        CaP.dfS = 0.0;


                                        CaP.fE += CaP.dfE;
                                        CaP.dfE = 0.0;


                                        #ifdef eglocal
                                        CA[ii*mx + ji].dfG = 0.0;
                                        #endif

                                        /*
                                        #ifdef eglocal
                                        CaPaux.Amin = 0.0;
                                        CaPaux.Amax = 0.0;
                                        #endif
                                        */

                                        if(fabs(CaP.state) >=2.0){
                                                CaP.Cd = CaP.Cddt;
                                                CaP.Cs = CaP.Csdt;
                                        }
                                        else{
                                                CaP.Cd = VfP.CLdt;
                                                CaP.Cs = kpart*VfP.CLdt;
                                        }


                                        CaPaux.state = 0;
                                        CaPaux.fS = 0.0;
                                        CaPaux.dfS = 0.0;
                                        CaPaux.theta = 0.0;
                                        CaPaux.L[0] = 0.0;
                                        CaPaux.L[1] = 0.0;
                                        CaPaux.L[2] = 0.0;
                                        CaPaux.L[3] = 0.0;
                                        CaPaux.dC[0] = 0.0;

                                        CaPaux.dC[1] = 0.0;
                                        CaPaux.Cd = VfP.CL;
                                        CaPaux.Cs = 0.0;
                                        CaPaux.Csdt = 0.0;




                                        CaPaux.fE = 0.0;
                                        CaPaux.dfE = 0.0;





                                }
                        }
                }
        }


}

void DiffusionLenght(){

                                velocity = DendriticVelocity();

                                #ifdef SeCACorrectionSphere
                                                velocity = 1.0/sqrt(2.0*MPI);
                                #endif

                                #ifdef GrowthDel
                                del_dl = SoluteDiffusionLength( fLInterpol(), velocity, NeInterpol() );
                                #endif
                                #ifdef ConvecDel
                                del_dl = SoluteDiffusionConvecLength( VfP.fL, fabs(VfP.VL - VfP.Vs), VfP.Ne );
                                #endif
                                #ifdef RanzMarshallDel
                                del_dl = RanzMarshallLenght(liquiddensity*fabs(clusters[ClusterP].V[1])*EquivalentD( clusters[ClusterP].size*dxCA*dyCA )/2.0/viscosity, viscosity/liquiddensity/DiffLiq, EquivalentD( clusters[ClusterP].size*dxCA*dyCA )/2.0);
                                #endif

                                //DelAdjustmentFunction();

                                #ifdef NewDeldl
                                double NePos = NeInterpol();
                                double fLPos = fLInterpol();
                                double Rf = pow(3.0*mmix*mmiy*dxCA*dyCA/4.0/M_PI/NePos, 1.0/3.0);
                                double Re = Rf*pow( 1.0 - fLPos , 1.0/3.0);
                                double Pe = velocity*Re/DiffLiq;
                                double Reynolds = liquiddensity*fabs(clusters[ClusterP].V[1])*Re/viscosity;
                                //double Reynolds = liquiddensity*fabs(Vref - VF[0].VL)*Re/viscosity;
                                double Schmidt = viscosity/liquiddensity/DiffLiq;

                                if(t > 0.0165)
                                                t=t;

                                del_dl = NewDelDL( Re, Rf, Pe, delC_Stag(Re, Rf, Reynolds,Schmidt) );
                                #endif

}

void DelAdjustmentFunction(){

                                DiffusionLenght();

                                /*
                                #ifdef DelAdjustment
                                        #ifndef SelectiveDelAdjustment
                                        if(clusters[ClusterP].fixed != 1){
                                                DelAdjustmentValue = 1.0;
                                                //DelAdjustmentValue = exp( - DelAdjustConstant*(fabs(BoundaryConditionHolder[ji].edVd) + fabs(clusters[ClusterP].V[1]) ) );
                                                DelAdjustmentValue = 0.1;
                                                if(DelAdjustmentValue > 1.0)
                                                        ii=ii;
                                                del_dl = del_dl*DelAdjustmentValue;

                                        }
                                        #endif
                                        #ifdef SelectiveDelAdjustment
                                        if(clusters[ClusterP].fixed != 1){

                                                ii=ii;
                                                if(ii != my - 1 && CaS.state == 0){

                                                        del_dl = DiffLiq/(velocity - (clusters[ClusterP].V[1] - BoundaryConditionHolder[ji].edVd));
                                                        if(DiffLiq/del_dl < fabs(velocity - clusters[ClusterP].V[1] + BoundaryConditionHolder[ji].edVd) ) {
                                                                del_dl = DiffLiq/(velocity - clusters[ClusterP].V[1] + BoundaryConditionHolder[ji].edVd);
                                                                ii=ii;
                                                        }



                                                }
                                                else if(ii != 0 && CaN.state == 0){
                                                        del_dl = DiffLiq/(velocity + (clusters[ClusterP].V[1] - BoundaryConditionHolder[ji].edVd));
                                                        //DelAdjustmentValue = exp(-DelAdjustConstant*(BoundaryConditionHolder[ji].edVd - clusters[ClusterP].V[1]));

                                                        if(DiffLiq/del_dl < (velocity + clusters[ClusterP].V[1] - BoundaryConditionHolder[ji].edVd) ){
                                                                del_dl = DiffLiq/(velocity + clusters[ClusterP].V[1] - BoundaryConditionHolder[ji].edVd);
                                                                ii=ii;
                                                        }


                                                }
                                                else{
                                                        del_dl = DiffLiq/(velocity);
                                                }
                                        }

                                        ii=ii;
                                        #endif
                                #endif
                                */

                                //del_dl = pow(3.0/4.0 , 1.0/3.0)*de*pow(Reynolds*Schmidt , -1.0/3.0)*0.87343322;
                                //del_dl = 1.0;


}

