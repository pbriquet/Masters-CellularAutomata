#include <math.h>
#include <stdio.h>
#include <stdlib.h>

//---------------------------------------------------------------------------
#pragma hdrstop

#include "PhaseTransf.h"
//---------------------------------------------------------------------------
#pragma package(smart_init)


#include "GlobalDefinitions.h"
#include "SetupDatabase.h"
#include "GlobalVariables.h"

#include "GrowthKnetics.h"
#include "InterpolationFunctions.h"
#include "CaMotion.h"
#include "DiffusionLenghts.h"
#include "ExtraFunctions.h"
#include "OutputFiles.h"








void dHDiffusive(){

        #ifndef ZeroDimensions
        double Tp,Tn,Ts;
        double Keqp, Keqn, Keqs;
        double Keqpn, Keqps;
        double fSp, fSn, fSs;

        // Cálculo de dHdiffusive

        double fluxn, fluxs;
        #endif


        for(ia=0;ia<mmacroy;ia++){
                for(ja=0;ja<mmacrox;ja++){

                        #ifdef VFHeatBalance
                        Tp = VfP.T;
                        fSp = VfP.fS;
                        TVF = Tp;
                        Keqp = fSp*Ksolid + (1.0 - fSp)*Kliquid;

                        if(ia!=0){

                                Tn = VfN.T;
                                fSn = VfN.fS;
                                TVF = Tn;
                                Keqn = (fSn)*Ksolid + (1.0 - fSn)*Kliquid;
                                Keqpn = 2.0*Keqn*Keqp/(Keqn + Keqp);

                                fluxn = Keqpn*(Tn - Tp)/dyVF;

                        }
                        if(ia!=mmacroy - 1){
                                Ts = VfS.T;
                                fSs = VfS.fS;
                                TVF = Ts;
                                Keqs = (fSs)*Ksolid + (1.0 - fSs)*Kliquid;
                                Keqps = 2.0*Keqs*Keqp/(Keqs + Keqp);

                                fluxs = Keqps*(Tp - Ts)/dyVF;
                        }
                        if(ia==0){
                                fluxn = -Hconvectionup*(Tboundary( Tp, Tfarfieldup, Keqp, Hconvectionup, dyVF) - Tfarfieldup );
                        }
                        if(ia==mmacroy - 1){
                                fluxs = Hconvectiondown*(Tboundary( Tp, Tfarfielddown, Keqp, Hconvectiondown, dyVF) - Tfarfielddown );
                        }

                        VfP.dH = (fluxn - fluxs)*dt/dyVF;
                        #endif

                        #ifdef ZeroDimensions
                                VfP.dH = CoolingRate*dt*soliddensity*Cp;
                        #endif

                        if( VfP.dH != 0.0 ){
                                ia=ia;
                        }

                        ia=ia;

                        #ifdef VirtualZeroDimensions
                                if( ia == mmay - 1){
                                                VfP.dH = CoolingRate*dt*soliddensity*Cp;
                                }
                        #endif

                        #ifdef HomogenousVFTemperature

                        #endif
                }
        }
}


void getCoefficientdfSMicro(){

                double fdPN;



                struct MicroMatrix cell;
                int k;

                for( k = 0; k < 1 + 2*dimmax ; k++){
                                afdCd[k] = 0.0;
                                afd[k] = 0.0;
                                vfd[k] = 0.0;
                                vfdCd[k] = 0.0;
                }

                fdP = 1.0 - CaP.fS;
                #ifdef Eutetic
                fdP -= CaP.fE;
                afdCdgamaE = -(Teut - Tf)/ml*dt/soliddensity;
                #endif
                vfdCd[0] =  fdP*CaP.Cd;
                afdCdgamaS = -kpart*CaP.Cd*dt/soliddensity;
                aTgamaS = dt/soliddensity/Cp*Lf;
                afdgamaS = -dt/soliddensity;
                afdCd[0] = 1.0;
                afd[0] = 1.0;
                vfd[0] = fdP;


                CdCA = CaP.Cd;
                CLVFDiffusion = VfP.CL;

                #ifdef CLInterpolforLiquidDiffusion
                        CLVFDiffusion = CLInterpol();
                #endif

                #ifdef CLequalCoforLiquidDiffusion
                        CLVFDiffusion = Co;
                #endif

                if(CdCA < VfP.CL)
                        ii = ii;

                for( k = 1 ; k <= 4 ; k ++ ){

                                if( k == 1 ){
                                                if(ii != 0)
                                                cell = CaN;
                                                else
                                                continue;
                                                //North
                                }
                                if( k == 2 ){
                                                if(ii != my - 1)
                                                cell = CaS;
                                                else
                                                continue;
                                                // South
                                }
                                if( k == 4 ){
                                                if(ji != 0)
                                                cell = CaW;
                                                else
                                                continue;
                                                // West
                                }
                                if( k == 3 ){
                                                if(ji != mx - 1)
                                                cell = CaE;
                                                else
                                                continue;
                                                // East
                                }

                                if( fabs(cell.state) >= 2){

                                                vfd[k] = 1.0 - cell.fS;
                                                #ifdef Eutetic
                                                vfd[k] -= cell.fE;
                                                #endif
                                                fdPN = 2.0/( 1.0/vfd[k] + 1.0/fdP );
                                                afdCd[k] = fdPN*InterliquidDiffusion*dt/dxCA/dxCA;
                                                vfdCd[k] = vfd[k]*cell.Cd;

                                                afdCd[0] -= afdCd[k]/fdP;
                                                afdCd[k] = afdCd[k]/vfd[k];

                                }


                }

                if(CA[ii*mx + ji].state >= 2){

                        if(fabs(VfP.fL - fLVFmin) > 1.0e-15){
                                #ifdef DiffusionThroughVFs
                                if(ii == ia*mmiy){
                                        #ifndef CLInterpolforLiquidDiffusion
                                        CLVFDiffusion = VfP.CL;
                                        #endif
                                        #ifdef BlockinterVFDiffusion
                                        CLVFDiffusion = CdCA;
                                        #endif
                                }
                                #endif

                                DelAdjustmentFunction();
                                if(del_dl <= 100.0)
                                                JDdlCA = dt*ExtraliquidDiffusion/del_dl;
                                else
                                                JDdlCA = 0.0;


                                #ifndef SeCACorrection
                                JDdlCA = 4.0*JDdlCA/dxCA;
                                #endif

                                #ifdef eglocal
                                #ifdef SeCACorrection
                                JDdlCA = JDdlCA*( 4.0*CaP.L[0]/(CaP.Amax - CaP.Amin) );
                                #endif
                                #endif

                                #ifdef SeCACorrectionDiagonal
                                JDdlCA = JDdlCA*sqrt(2.0);
                                #endif

                                #ifndef DiffusionThroughVFs
                                JDdlCA = JDdlCA*(CdCA - CLVFDiffusion);
                                #endif

                                #ifdef DiffusionThroughVFs
                                if(ii != ia*mmiy)
                                JDdlCA = JDdlCA*(CdCA - CLVFDiffusion);
                                #endif



                                #ifdef CdlEqualFluxinDiffusion




                                if(del_dl <= 100.0){
                                Rdl = ExtraliquidDiffusion/del_dl;
                                Rp = InterliquidDiffusion/dxCA*2.0*fdP;
                                CdlCA = (Rp*CdCA + Rdl*CLVFDiffusion)/(Rp + Rdl);
                                }
                                else
                                                CdlCA = CLVFDiffusion;
                                #ifdef CdlCAConstante
                                CdlCA = CdlCAConst;
                                #endif
                                
                                JDdlCA = 4.0*Rdl*dt/dyCA*(CdlCA - CLVFDiffusion);



                                if(ia == 0)
                                        ia=ia;
                                if(ia == 1)
                                        ia=ia;

                                #ifdef KeepClinCA
                                CaP.Cl = CdlCA;
                                #endif
                                #endif

                                #ifdef NSidesInLiquidDiffusion
                                int Sides = 0.0;
                                if( ii != 0.0 && CaN.state == 0)
                                                Sides ++;
                                if( ii != my - 1 && CaS.state == 0)
                                                Sides ++;
                                if( ji != mx - 1 && CaE.state == 0)
                                                Sides ++;
                                if( ji != 0.0 && CaW.state == 0)
                                                Sides ++;

                                
                                JDdlCA = JDdlCA*Sides/4.0;

                                if(JDdlCA < 0)
                                        ii=ii;

                                #endif

                                #ifdef STEREOMICRO
                                JDdlCA = (4.0/M_PI)*JDdlCA;
                                #endif

                                if(JDdlCA < 0.0)
                                                ii=ii;
                        }
                        else{
                                JDdlCA = 0.0;
                        }

                }



                #ifdef inWallThereisInterExtraDiffusion

                #endif

                JDdlVF += JDdlCA/mmix/mmiy;

                SUMafdvfd = 0.0;
                SUMafdCdvfdCd = 0.0;

                for( int k = 0; k < 1 + 2*dimCA ; k++){
                                SUMafdvfd += afd[k]*vfd[k];
                                SUMafdCdvfdCd += afdCd[k]*vfdCd[k];
                }

                SUMafdCdvfdCd -= JDdlCA;

                if(afd[0] <= 0.0 || afdCd[0] < 0.0){

                        ii=ii;
                        printf("/n/nAlert! afdCd[0] < 0!");
                        //system("pause");
                }

                if(afdCd[0] - JDdlCA/(CdCA - CLVFDiffusion + 0.0000001)/fdP < 0.0)
                                ii=ii;




}


void getCoefficientsdCLMacro(){



        SUMaflvfl = 0.0;
        SUMaflClvflCl = 0.0;

        for( int k = 0; k < 1 + 2*dimmax ; k++){
                                aflCl[k] = 0.0;
                                afl[k] = 0.0;
                                vfl[k] = 0.0;
                                vflCl[k] = 0.0;
        }

        struct MacroMatrix cell;

        int k = 0;
        cell = VfP;
        //vflCl[k] = cell.fL*cell.CL;
        //vfl[k] = cell.fL;
        vflCl[k] = cell.fLCA*cell.CL;
        vfl[k] = cell.fLCA;
        afl[k] = 1.0;
        aflCl[k] = 1.0;

        if(ia != 0){
                cell = VfN;
                k = 1;

                vflCl[k] = cell.fL*cell.CL;
                vfl[k] = cell.fL;

                #ifndef CLVFDiffusionDiffusion
                aflCl[k] = 0.0;
                afl[k] = 0.0;
                #endif

                #ifdef CLVFDiffusionDiffusion
                aflCl[k] = 2.0*VfP.fL/(cell.fL + VfP.fL)*DiffLiq*dt/dxVF/dxVF;
                afl[k] = 0.0;

                #ifdef DiffusionDistanceColumnar
                                                aflCl[k] = aflCl[k]/((VfP.fL + 1.0)/2.0);
                                #endif

                aflCl[0] -= aflCl[k]/VfP.fL*cell.fL;
                afl[0] -= afl[k];
                #endif

        }
        if( ia != mmay - 1){
                cell = VfS;
                k = 2;

                vflCl[k] = cell.fL*cell.CL;
                vfl[k] = cell.fL;

                #ifndef CLVFDiffusionDiffusion
                aflCl[k] = 0.0;
                afl[k] = 0.0;
                #endif

                #ifdef CLVFDiffusionDiffusion

                aflCl[k] = 2.0*VfP.fL/(cell.fL + VfP.fL)*DiffLiq*dt/dxVF/dxVF;
                afl[k] = 0.0;

                                #ifdef DiffusionDistanceColumnar
                                                aflCl[k] = aflCl[k]/((VfP.fL + 1.0)/2.0);
                                #endif

                aflCl[0] -= aflCl[k]/VfP.fL*cell.fL;
                afl[0] -= afl[k];
                #endif


        }

        if(ja != 0){

                cell = VfW;
                k = 4;

                vflCl[k] = cell.fL*cell.CL;
                vfl[k] = cell.fL;

                #ifndef CLVFDiffusionDiffusion
                aflCl[k] = 0.0;
                afl[k] = 0.0;
                #endif

                #ifdef CLVFDiffusionDiffusion
                aflCl[k] = 2.0*VfP.fL/(cell.fL + VfP.fL)*DiffLiq*dt/dxVF/dxVF;
                afl[k] = 0.0;

                aflCl[0] -= aflCl[k]/VfP.fL*cell.fL;
                afl[0] -= afl[k];
                #endif


        }
        if(ja != mmax - 1){
                cell = VfE;
                k = 3;

                vflCl[k] = cell.fL*cell.CL;
                vfl[k] = cell.fL;

                #ifndef CLVFDiffusionDiffusion
                aflCl[k] = 0.0;
                afl[k] = 0.0;
                #endif

                #ifdef CLVFDiffusionDiffusion
                aflCl[k] = 2.0*VfP.fL/(cell.fL + VfP.fL)*DiffLiq*dt/dxVF/dxVF;
                afl[k] = 0.0;

                aflCl[0] -= aflCl[k]/VfP.fL*cell.fL;
                afl[0] -= afl[k];

                #endif

        }


        for( int k = 0; k < 1 + 2*dimVF ; k++){
                SUMaflvfl += afl[k]*vfl[k];
                SUMaflClvflCl += aflCl[k]*vflCl[k];
        }

        if(afl[0] <= 0.0 || aflCl[0] < 0.0){

                        ii=ii;
                        printf("/n/nAlert! aflCl[0] < 0!");
                        system("pause");
        }


}


void dfSMacro(){

                        VF[ia*mmax + ja].dfS = VF[ia*mmax + ja].dfS/mmix/mmiy;

                        #ifdef Eutetic
                        VF[ia*mmax + ja].dfE = VF[ia*mmax + ja].dfE/mmix/mmiy;
                        #endif

}

void dTMacro(){



        VfP.dT = VfP.dH/soliddensity/Cp + Lf/Cp*VfP.dfS;
        #ifdef Eutetic
        VfP.dT += Lf/Cp*VfP.dfE;
        #endif

        if( VfP.dT >= 0.0 )
                ia=ia;

        if( VfP.dT <= 0.0 )
                ia=ia;

        #ifdef HomogenousVFTemperature
        if(ia != mmay - 1){
                VfP.dT = VF[mmay - 1].dT;
        }
        #endif


}

void dCLMacro(){

double fLdt;
        VfP.fL = VfP.fLCA;
        getCoefficientsdCLMacro();
        if(JDdlVF != 0.0){
                ia=ia;
        }
        if(JldVl != 0.0)
                ii=ii;

        CdJldVl = CdJldVl*dt*dxCA/dyVF/dxVF;
        CdJldVs = CdJldVs*dt*dxCA/dyVF/dxVF;
        JldVl = JldVl*dt*dxCA/dyVF/dxVF;
        JldVs = JldVs*dt*dxCA/dyVF/dxVF;


        //if( fabs(JldVs) < fabs(0.001*dxCA/dyVF*VfP.Vs))

        if( JldVl  > 0.0){
                /*
                if( ia == mmay - 1 || ia != mmay - 1 && VfS.fL == fLVFmin){
                        if(JldVl > 0.0){
                                JldVl = 0.0;
                                JldVs = 0.0;
                        }
                }
                */
                /*

                if(fabs(JldVl - JldVs) < (1.0 - VfP.fL - VfP.fS - VfP.fE)*fabs(VfP.Vd - VfP.Vs)*dt*fabs(AldOpenDown + AldOpenUp)*dxCA/dxVF/dyVF){
                        JldVl = 0.0;
                        JldVs = 0.0;
                }
                */


                if( fabs(AldOpenDown + AldOpenUp) < 1.0 ){  //fabs(AldOpen)  < 0.1*dxCA/dxVF &&
                        if(ia != mmay-1){
                        if(JldVl != elvlyn[ia+1]*dt/dyVF){
                                JldVl = elvlyn[ia+1]*dt/dyVF;
                        }
                        }
                        if(ia == mmay - 1){
                                JldVl = 0.0;
                                JldVs = 0.0;
                        }
                        /*
                        if(ia != mmay - 1 && elvlyn[ia+1] == 0.0){
                                JldVl = 0.0;
                                JldVs = 0.0;
                        }
                        */

                }
                else
                        ii=ii;





        }





        double StereoFactor;

        #ifdef STEREOMACROGROWTH
        CdJdlVF = CdJdlVF*(4.0/M_PI);
        JdlVF = JdlVF*(4.0/M_PI);
        #endif


        #ifndef STEREOMICRO
        #ifdef STEREOMACRO
        JDdlVF = JDdlVF*(4.0/M_PI); // Stereological Transformation
        #endif
        #endif


        #ifdef MotionAlgorithm
        #ifdef STEREOOneSphereMACRO
        StereoFactor = M_PI/4.0*fabs(AldOpen);
        if(StereoFactor >= 1.0e-10){
                CdJldVl = CdJldVl*StereoFactor;
                CdJldVs = CdJldVs*StereoFactor;
                JldVs = JldVs*StereoFactor;
                JldVl = JldVl*StereoFactor;
        }
        if(StereoFactor >= M_PI/8.0 - 0.001)
                ii=ii;


        #endif
        #ifdef STEREOOneOctahedronMACRO
        StereoFactor = fabs(AldOpen);
        if(StereoFactor >= 1.0e-10){
                CdJldVl = CdJldVl*StereoFactor;
                CdJldVs = CdJldVs*StereoFactor;
                JldVs = JldVs*StereoFactor;
                JldVl = JldVl*StereoFactor;
        }

        #endif


        elvlyn[ia] -= JldVl*dyVF/dt;

        if(ia != mmay - 1)
                elvlyn[ia] += elvlyn[ia + 1];

        if(ia == 0)
                elvlyn[ia] = 0.0;

        if(ia != mmay - 1 && VfP.fL == 1.0 && elvlyn[ia-1] != 0.0)
                elvlyn[ia] = 0.0;


        if(elvlyn[ia] < 0.0){
                //PrintSurfaceCd();
                ii=ii;
                //elvlyn[ia] = 0.0;
                /*
                for(RunInverseCAinVFLine){
                        for(RunCAinVFCol){
                                if(CaP.state != 0)
                                        dfSMicro();
                        }
                }
                */
                /*
                FILE * Errorelvl;
                Errorelvl = fopen("Erelvl.dat","w");
                fprintf(Errorelvl, "\nError! t = %f", t);
                */
                /*
                for(RunCAinVFLine){
                for(RunCAinVFCol){
                        if(CaP.state != 0.0)
                                fprintf(Errorelvl, "%.20f;", clusters[ClusterP].V[1]);
                        else
                                fprintf(Errorelvl, "%.20f;", 0.0);

                }
                fprintf(Errorelvl, "\n");
                }

                fprintf(Errorelvl, "\n\n\n");
                for(RunCAinVFLine){
                for(RunCAinVFCol){
                        fprintf(Errorelvl, "%3.d;", CaP.state);

                }
                fprintf(Errorelvl, "\n");
                }

                fprintf(Errorelvl, "\n\n\n");
                for(RunCAinVFLine){
                for(RunCAinVFCol){
                        fprintf(Errorelvl, "%3.d;", clusters[ (int) fabs(CaP.state) - 2].cluster);

                }
                fprintf(Errorelvl, "\n");
                }

                PrintGrains();
                */
        }
        
        if(elvlyn[ia]!=0)
                ii=ii;

        #ifdef CLConvectiveTermCalculation
        double CLConvectiveTerm;
        double fLyn, CLyn, fLys, CLys;
        double vLyn, vLys;

        #ifdef fLCLupwind
        if(ia != 0){
                if(elvlyn[ia] >= 0.0){
                        fLyn = VfP.fL;
                        CLyn = VfP.CL;
                }
                if(elvlyn[ia] < 0.0){
                        fLyn = VfN.fL;
                        CLyn = VfN.CL;
                }
        }
        else{
                fLyn = 0.0;
                CLyn = 0.0;
        }
        if(ia != mmay - 1){
                if(elvlyn[ia + 1] >= 0.0){
                        fLys = VfS.fL;
                        CLys = VfS.CL;
                }
                if(elvlyn[ia+1] < 0.0){
                        fLys = VfP.fL;
                        CLys = VfP.CL;
                }
        }
        else if(ia == mmay - 1){
                fLys = 0.0;
                CLys = 0.0;
        }



        else{
                fLys = 0.0;
                CLys = 0.0;
        }
        #endif

        /*
        #ifdef vLupwind
        vLyn = ia!=0 ? VfP.VL : 0.0;
        if(ia != mmay - 1)
                vLys = VfS.VL;
        else
                vLys = 0.0;
        #endif

        #ifdef vLCentral
        vLyn = ia!=0 ? (VfN.VL + VfP.VL)/2.0 : 0.0 ;
        if(ia != mmay - 1)
                vLys = (VfS.VL + VfP.VL)/2.0;
        else
                vLys = 0.0;
        #endif
        */

        if(ia != mmay - 1){
                vLys = elvlyn[ia+1];
        }
        else{
                vLys = 0.0;
        }
        vLyn = elvlyn[ia];

        #ifdef ReferenceTest
        vLys = -Vref;
        vLyn = -Vref - JldVl*dyVF/dt;
        CLyn = VF[0].CL;
        CLys = VF[0].CL;
        #endif


        #ifdef CLdtConservative
        CLConvectiveTerm = -(vLyn*CLyn - vLys*CLys)*dt/dyVF;
        #endif
        #ifdef CLdtSemiConservative
        CLConvectiveTerm = -(vLyn*CLyn - vLys*CLys)*dt/dyVF;
        #endif
        #ifdef CLdtNonConservative
        CLConvectiveTerm = -vLyn*fLyn*(CLyn - CLys)*dt/dyVF;
        #endif



        #ifdef CdJldExtraTerm
                                if(CdJldVl - CdJldVs != 0.0)
                                        ii=ii;

                                if(CdJldVs != 0.0)
                                        ii=ii;

                #ifdef JldVsDiscrete
                
                #endif
                CdJdlVF -= -CdJldVl + CdJldVs;
                JdlVF -= JldVs;
        #endif

        if(CLConvectiveTerm != 0.0)
                ii=ii;

        SUMaflClvflCl += CLConvectiveTerm;

        #endif

        // This Refresh Value of fL can become an error
        #ifdef RefreshfLdtvalueConvective
        //SUMaflvfl += -(vLyn*fLyn - vLys*fLys)*dt/dyVF;

        fLdt = SUMaflvfl - JdlVF;

        #endif
        #endif

        /*
        #ifdef TruncarfLVelocidade
        if( fLdt < fLVFmin){
                fLdt = fLVFmin;
                if(JdlVF != 0.0){
                                CdJdlVF = -CdJdlVF/JdlVF*(fLdt - VfP.fL);
                                JdlVF = fLdt - VfP.fL;
                }
        }
        #endif
        */

        #ifdef TruncarfLdt

        if( fLdt < fLVFmin){
                fLdt = fLVFmin;

                if(JdlVF == 0.0){
                                JDdlVF = 0.0;
                                CdJdlVF = 0.0;
                }

                if(JdlVF != 0.0){
                                CdJdlVF = -CdJdlVF/JdlVF*(fLdt - VfP.fL);
                                JDdlVF = -JDdlVF/JdlVF*(fLdt - VfP.fL);
                                JdlVF = -(fLdt - VfP.fL);
                }

                /*
                JDdlVF = 0.0;
                CdJdlVF = 0.0;
                JdlVF = 0.0;
                */

        }
        if(fabs(VfP.fL-fLVFmin) < 1.0e-15){
                JDdlVF = 0.0;
                CdJdlVF = 0.0;
                JdlVF = 0.0;
        }
        /*
        if( fLdt < fLVFmin ){
                fLdt = fLVFmin;
        }
        */
        #endif


        #ifdef CLdtConservative
        if(VfP.CLdt - VfP.CL != 0.0)
                ia=ia;
        VfP.CLdt = (SUMaflClvflCl - CdJdlVF + JDdlVF)/(fLdt);
        #endif
        #ifdef CLdtSemiConservative
        VfP.CLdt = (SUMaflClvflCl - CdJdlVF + JDdlVF + (JdlVF)*VfP.CL)/(VfP.fL);
        #endif
        /*
        #ifdef CLdtNonConservative
        VfP.CLdt = (SUMaflClvflCl - (CdJdlVF - VfP.CL*JdlVF) + JDdlVF)/(VfP.fL);
        #endif
        */
        if(JdlVF != 0)
                ia=ia;

        if(JDdlVF != 0)
                ia=ia;
        /*
        if( fLdt <= fLVFmin){
                fLdt = fLVFmin;
                VfP.CLdt = VfP.CL + VfP.fL/fLVFmin*(VfP.CLdt - VfP.CL);
        }
        */
        if( fLdt > fLVFmin){

        }

        fLdt = VfP.fLCA;

        VfP.dfL = fLdt - VfP.fL;


        if( VfP.CLdt < 0.0 ){
                ia=ia;
        }

        if( VfP.CLdt != Co ){
                ia=ia;
        }

        if( VF[0].CLdt > VF[1].CL){
                ia=ia;
        }
        if( VfP.CLdt < Co - 0.01 ){
                ia=ia;
                //printf(" CL(%d) ", ia);
        }

        if(VfP.CLdt < Co)
                ia=ia;

        if( VfP.CLdt > 20.0)
                ia=ia;

        if( VfP.CLdt < VfP.CL )
                ia=ia;


        if( VF[0].CLdt < VF[0].CL )
                ia=ia; 
        //VfP.CLdt = Co;
        // Condição do Gandin


}


void dfLMacro(){

#ifndef DeterministicModel

        #ifndef eglocal

        #ifdef JldVsDiscrete
        if(fabs(CAaux[ii*mx+ji].state) >= 2 && CaP.state == 0.0){
                CdJldVs += CAaux[ii*mx+ji].Cl*dyCA*dxCA/dyVF/dxVF;
        }
        if(CaP.state != 0 && CAaux[ii*mx+ji].state == 0 && fabs(CA[ii*mx + ji].theta) <= M_PI/4.0 + dtheta){
                CdJldVs -= CA[ii*mx+ji].Cl*dyCA*dxCA/dyVF/dxVF;
        }
        #endif

        if(fabs(CA[ii*mx + ji].state) >= 2 && fabs(CA[ii*mx + ji].theta) > M_PI/4.0 + dtheta){



                CA[ii*mx + ji].theta -= 10.0;
                JdlCA = 1.0;
                JdlVF += JdlCA/mmix/mmiy;
                #ifndef CdlEqualFluxinDiffusion
                CdJdlVF += Cld*JdlCA/mmix/mmiy;
                #endif
                #ifdef CdlEqualFluxinDiffusion
                CdJdlVF += CaP.Cd*JdlCA/mmix/mmiy;
                #endif

                clusters[ ClusterP].fS = clusters[ ClusterP ].fS*clusters[ ClusterP ].size/(clusters[ClusterP].size + 1);
                clusters[ ClusterP].size ++;

                if(clusters[ ClusterP ].size > mx*my)
                        ii=ii;



        }


        #endif

        #ifdef eglocal
        if(fabs(CA[ii*mx + ji].state) >= 2 && fabs(CA[ii*mx + ji].theta) > M_PI/4.0 + dtheta){

                CA[ii*mx + ji].theta -= 10.0;
                CalculateAmax();
                clusters[ clusters[(int) fabs(CA[ii*mx + ji].state) - 2].cluster].size ++;

        }

        if(CA[ii*mx + ji].state >= 2){
                fGCAdt = (2.0*CaP.L[0]*CaP.L[0] - CaP.Amin)/(CaP.Amax - CaP.Amin);
                fGCA = fGCAdt - CA[ii*mx + ji].dfG;
                if(CA[ii*mx + ji].dfG != 0.0)
                                ii=ii;
                JdlCA = CA[ii*mx + ji].dfG;

                #ifdef NSidesInGrainGrowth
                                int Sides = 0.0;
                                if( ii != 0.0 && CaN.state == 0)
                                                Sides ++;
                                if( ii != my - 1 && CaS.state == 0)
                                                Sides ++;
                                if( ji != mx - 1 && CaE.state == 0)
                                                Sides ++;
                                if( ji != 0.0 && CaW.state == 0)
                                                Sides ++;

                                JdlCA = JdlCA*Sides/4.0;
                #endif

                if(JdlCA != 0.0)
                ii =ii;

                #ifdef CdlEqualFluxinDiffusion
                fdP = 1.0 - CaP.fS - CaP.fE;
                velocity = DendriticVelocity();
                CLVFDiffusion = VfP.CL;
                CdCA = CaP.Cd;
                #ifdef SeCACorrectionSphere
                velocity = 1.0/sqrt(2.0*MPI);
                #endif
                DiffusionLength();
                Rdl = ExtraliquidDiffusion/del_dl;
                Rp = InterliquidDiffusion/dxCA*2.0*fdP;
                CdlCA = (Rp*CdCA + Rdl*CLVFDiffusion)/(Rp + Rdl);

                #endif


                JdlVF += JdlCA/mmix/mmiy;
                CdJdlVF += Cld*JdlCA/mmix/mmiy;

        }


        #endif

#endif


}


#ifdef Eutetic
void dfSMicro(){
                double TCAt;
                double TCAtdt;
                #ifndef ScheilMicrofS
                double gamaS;
                double fSCAtdt;

                double gamaE;
                double fECAtdt;
                double fdCA;
                double fdCAtdt;
                double gamaSgamaE;


                double vs = 0.0;

                double Tp, Tn, Ts, TConvecTerm;
                double Cdp, Cds, Cdn, CdConvecTerm;

                double a,b,c;
                double gamaS1;
                double gamaS2;

                if(t >= 0.1965 && t <= 0.1965 + dt)
                                ii==ii;  // ***



                
                if( fabs(CaP.state) >= 2 ){





                dHCA = dHInterpol();



                #ifdef TInterpolinPhaseTransformation
                TCA = TInterpol();
                #endif

                #ifdef TinEquilibriuminPhaseTransformation
                if(CaP.fS != 0.0){
                                TCA = Tf + CaP.Cd*ml;
                }
                else{
                                TCA = TInterpol();
                }
                #endif

                CdCA = CA[ii*mx + ji].Cd;
                fSCA = CA[ii*mx + ji].fS;
                fECA = CA[ii*mx + ji].fE;
                fdCA = 1.0 - fSCA - fECA;

                #ifdef eglocal
                                fGCAdt = (2.0*CaP.L[0]*CaP.L[0]- CaP.Amin)/(CaP.Amax - CaP.Amin);
                                if(fGCAdt < 0.0)
                                                fGCAdt = 0.0;
                                else if(fGCAdt > 1.0)
                                                fGCAdt = 1.0;
                                fGCA = fGCAdt - CA[ii*mx + ji].dfG;
                #endif

                #ifdef MotionAlgorithm
                                #ifndef MartoranoLiquidBoundaryCondition
                                vs = clusters[ClusterP].V[1];
                                #endif

                                #ifndef HoldBoundaryCondition
                                vs = -(1.0 - VfP.fS - VfP.fE - VfP.fL)*VfP.Vd + (1.0 - CaP.fS - CaP.fE)*clusters[ClusterP].V[1];
                                //vs = 0.001 - 0.005;
                                #endif

                                #ifdef HoldBoundaryCondition

                                if(BoundaryConditionHolder[ji].on == 0){
                                        BoundaryConditionHolder[ji].on = 1;
                                        BoundaryConditionHolder[ji].vcluster = clusters[ClusterP].V[1];



                                        if(CaP.fS + CaP.fE >= fSmax){
                                                BoundaryConditionHolder[ji].edVd = 0.0;
                                        }
                                        else if(clusters[ClusterP].fixed == 0){
                                                /*
                                                if(VfP.fL == 1.0)
                                                        BoundaryConditionHolder[ji].edVd = 0.0;
                                                else
                                                BoundaryConditionHolder[ji].edVd = (1.0 - VfP.fS - VfP.fE - VfP.fL)*(VfP.Vd - VfP.Vs)/(1.0 - VfP.fL);
                                                */

                                                if(fLInterpol() != 1.0){
                                                        BoundaryConditionHolder[ji].edVd = edVdVsInterpol()/(1.0 - fLInterpol());
                                                        //BoundaryConditionHolder[ji].edVd = (1.0 - VfP.fS - VfP.fE - VfP.fL)*(VfP.Vd - VfP.Vs)/(1.0 - VfP.fL);
                                                        //BoundaryConditionHolder[ji].edVd = 0.0;
                                                }
                                                else{
                                                        if( edVdVsInterpol() == 0.0 )
                                                                BoundaryConditionHolder[ji].edVd = 0.0;
                                                        else
                                                                BoundaryConditionHolder[ji].edVd = edVdVsInterpol()/(dxCA*dyCA/dxVF/dyVF);
                                                }
                                           
                                        }

                                        if( clusters[ClusterP].fixed == 1 || BoundaryConditionHolder[ji].vcluster == 0.0)
                                                BoundaryConditionHolder[ji].edVd = 0.0;

                                        //BoundaryConditionHolder[ji].edVd = 0.001 - 0.005;
                                        //vs = -BoundaryConditionHolder[ji].edVd;

                                        //if(BoundaryConditionHolder[ji].edVd != 0.0){
                                                AldOpenDown += 1.0;
                                                AldOpen += 1.0;
                                        //}
                                }


                                #endif
                #endif

                getCoefficientdfSMicro();

                if(BoundaryConditionHolder[ji].edVd < 0.0)
                        ii=ii;
                #ifdef MotionAlgorithm
                                  // It will be the (- Boundary Condition)

                                // edvd Boundary Condition Options

                                // Alternatives are:
                                // 1 -  #ifndef MartoranoLiquidBoundaryCondition
                                //      It means that edvd is not calculated by the VF Kv, and boundary condition is just -Vs (inverted signal for algorithm)
                                // 2 - #ifndef HoldBoundaryCondition
                                //      Constant is calculated in each CA, taking information from the VF in which CA belongs
                                // 3 - #ifdef HoldBoundaryCondition
                                //      Will Hold the Boundary Condition from interfface below through the CA Grain Column



                                // End of edvd Boundary Condition Options

                #ifdef TEqualFluxinDiffusion
                                // Convection Term is grouped with Tp in TCA. Tp will hold value at center for debug porpuses.

                                Tp = TCA;

                                if(BoundaryConditionHolder[ji].edVd != 0.0){
                                                if(ii!=my-1){
                                                                if(fabs(CaS.state) >= 2){
                                                                                //Ts = Tf + ml*CaS.Cd; // Up-Wind
                                                                                Ts = TInterpoly(0); // Up-Wind
                                                                }
                                                                else{
                                                                                /*
                                                                                DiffusionLenght();
                                                                                DelAdjustmentFunction();
                                                                                if(del_dl <= 100.0){
                                                                                                Rdl = ExtraliquidDiffusion/del_dl;
                                                                                                Rp = InterliquidDiffusion/dxCA*2.0*fdP;
                                                                                                CdlCA = (Rp*CdCA + Rdl*CLVFDiffusion)/(Rp + Rdl);
                                                                                }
                                                                                else
                                                                                                CdlCA = CLVFDiffusion;

                                                                                #ifdef CdlCAConstante
                                                                                CdlCA = CdlCAConst;
                                                                                #endif
                                                                                */
                                                                                //Ts = Tf + ml*CdlCA; // This works because CdlCA is always calculated in boundaries
                                                                                Ts = TInterpoly(0);
                                                                }
                                                }
                                                else{
                                                                //Ts = Tp;  // Didn't know what to put in this value
                                                                Ts = TInterpoly(0);
                                                }

                                                if(ii != 0){
                                                                if(fabs(CaN.state) >= 2){
                                                                                //Tn = Tf + ml*CaP.Cd; // Up-wind
                                                                                Tn = TInterpoly(1); // Interpolation
                                                                }
                                                                else{

                                                                                /*
                                                                                //DiffusionLength();
                                                                                DelAdjustmentFunction();
                                                                                if(del_dl <= 100.0){
                                                                                                Rdl = ExtraliquidDiffusion/del_dl;
                                                                                                Rp = InterliquidDiffusion/dxCA*2.0*fdP;
                                                                                                CdlCA = (Rp*CdCA + Rdl*CLVFDiffusion)/(Rp + Rdl);
                                                                                }
                                                                                else
                                                                                                CdlCA = CLVFDiffusion;

                                                                                #ifdef CdlCAConstante
                                                                                CdlCA = CdlCAConst;
                                                                                #endif
                                                                                */
                                                                                //Tn = Tf + ml*CdlCA;
                                                                                Tn = TInterpoly(1);
                                                                }
                                                }

                                                else{
                                                                //Tn = Tp; // Didn't know what to put in this value
                                                                Tn = TInterpoly(1);
                                                }
                                                TConvecTerm = -BoundaryConditionHolder[ji].edVd*dt/dyCA*(Tn-Ts);
                                                if(TConvecTerm != 0.0)
                                                                ii=ii;

                                                TCA += TConvecTerm;

                                }
                                else{
                                                TConvecTerm = 0.0;
                                }
                #endif


                                Cdp = CdCA;

                                //if(BoundaryConditionHolder[ji].edVd != 0.0){
                                                if(ii!=my-1){
                                                                if(fabs(CaS.state) >= 2){
                                                                                Cds = CaS.Cd; // Up-Wind
                                                                }
                                                                else{
                                                                                //DiffusionLength();
                                                                                DelAdjustmentFunction();
                                                                                #ifdef CdlEqualFluxinDiffusion
                                                                                if(del_dl <= 100.0){
                                                                                                Rdl = ExtraliquidDiffusion/del_dl;
                                                                                                Rp = InterliquidDiffusion/dxCA*2.0*fdP;
                                                                                                CdlCA = (Rp*CdCA + Rdl*CLVFDiffusion)/(Rp + Rdl);
                                                                                }
                                                                                else
                                                                                                CdlCA = CLVFDiffusion;
                                                                                #endif

                                                                                #ifdef CdlCAConstante
                                                                                CdlCA = CdlCAConst;
                                                                                #endif

                                                                                Cds = Cdl; // This works because CdlCA is always calculated in boundaries
                                                                                #ifdef CdJldExtraTerm

                                                                                CdJldVl += Cds*(BoundaryConditionHolder[ji].edVd + BoundaryConditionHolder[ji].vcluster );
                                                                                #ifndef JldVsDiscrete
                                                                                CdJldVs += Cds*( BoundaryConditionHolder[ji].vcluster );
                                                                                JldVs += ( BoundaryConditionHolder[ji].vcluster  );
                                                                                #endif
                                                                                JldVl += (BoundaryConditionHolder[ji].edVd +  BoundaryConditionHolder[ji].vcluster );


                                                                                // Sum for elvlyn
                                                                                /*
                                                                                #ifndef STEREOOneSphereMICRO
                                                                                elvlyn[ia*mmax + ja] -= dxCA*(-vs + clusters[ClusterP].V[1])/dxVF;
                                                                                #endif
                                                                                #ifdef STEREOOneSphereMICRO
                                                                                elvlyn[ia*mmax + ja] -= dxCA*(-vs + clusters[ClusterP].V[1])/dxVF*(M_PI/8.0);
                                                                                #endif
                                                                                */

                                                                                #endif


                                                                }
                                                }
                                                else{

                                                                //DiffusionLength();
                                                                DelAdjustmentFunction();

                                                                #ifdef CdlEqualFluxinDiffusion
                                                                if(del_dl <= 100.0){
                                                                                                Rdl = ExtraliquidDiffusion/del_dl;
                                                                                                Rp = InterliquidDiffusion/dxCA*2.0*fdP;
                                                                                                CdlCA = (Rp*CdCA + Rdl*CLVFDiffusion)/(Rp + Rdl);
                                                                                }
                                                                                else
                                                                                                CdlCA = CLVFDiffusion;
                                                                #endif

                                                                #ifdef CdlCAConstante
                                                                CdlCA = CdlCAConst;
                                                                #endif

                                                                Cds = CaP.Cd; // This works because CdlCA is always calculated in boundaries

                                                                #ifndef ReferenceTest

                                                                #ifdef CdJldExtraTerm
                                                                CdJldVl += Cds*(BoundaryConditionHolder[ji].edVd + BoundaryConditionHolder[ji].vcluster );
                                                                #ifndef JldVsDiscrete
                                                                CdJldVs += Cds*( BoundaryConditionHolder[ji].vcluster );
                                                                JldVs += ( BoundaryConditionHolder[ji].vcluster  );
                                                                #endif
                                                                JldVl += (BoundaryConditionHolder[ji].edVd +  BoundaryConditionHolder[ji].vcluster );


                                                                #endif

                                                                // Sum for elvlyn
                                                                /*
                                                                #ifndef STEREOOneSphereMICRO
                                                                elvlyn[ia*mmax + ja] -= dxCA*(-vs + clusters[ClusterP].V[1])/dxVF;
                                                                #endif
                                                                #ifdef STEREOOneSphereMICRO
                                                                elvlyn[ia*mmax + ja] -= dxCA*(-vs + clusters[ClusterP].V[1])/dxVF*(M_PI/8.0);
                                                                #endif
                                                                */
                                                                #endif
                                                }
                                                //if(ii == (ia + 1)*mmiy - 1){
                                                        //CdJldExtra += Cdp*vs*dxCA*dt/dxVF/dyVF;
                                                //}
                                                if(ii != 0){
                                                                if(fabs(CaN.state) >= 2){
                                                                                Cdn = CaP.Cd; // Up-wind
                                                                }
                                                                else{
                                                                                //DiffusionLength();
                                                                                DelAdjustmentFunction();
                                                                                #ifdef CdlEqualFluxinDiffusion
                                                                                if(del_dl <= 100.0){
                                                                                                Rdl = ExtraliquidDiffusion/del_dl;
                                                                                                Rp = InterliquidDiffusion/dxCA*2.0*fdP;
                                                                                                CdlCA = (Rp*CdCA + Rdl*CLVFDiffusion)/(Rp + Rdl);
                                                                                }
                                                                                else
                                                                                                CdlCA = CLVFDiffusion;
                                                                                #endif

                                                                                #ifdef CdlCAConstante
                                                                                CdlCA = CdlCAConst;
                                                                                #endif

                                                                                Cdn = Cdl;

                                                                                #ifdef CdJldExtraTerm
                                                                                CdJldVl += -Cdn*(BoundaryConditionHolder[ji].edVd + BoundaryConditionHolder[ji].vcluster);
                                                                                #ifndef JldVsDiscrete
                                                                                CdJldVs += -Cdn*( BoundaryConditionHolder[ji].vcluster);
                                                                                JldVs += -( BoundaryConditionHolder[ji].vcluster);
                                                                                #endif
                                                                                JldVl += -(BoundaryConditionHolder[ji].edVd + BoundaryConditionHolder[ji].vcluster);


                                                                                /*
                                                                                #ifndef STEREOOneSphereMICRO
                                                                                elvlyn[ia*mmax + ja] += dxCA*(-vs + clusters[ClusterP].V[1])/dxVF;
                                                                                #endif
                                                                                #ifdef STEREOOneSphereMICRO
                                                                                elvlyn[ia*mmax + ja] += dxCA*(-vs + clusters[ClusterP].V[1])/dxVF*(M_PI/8.0);
                                                                                #endif
                                                                                */

                                                                                #endif

                                                                }
                                                }
                                                else{
                                                                //DiffusionLength();
                                                                DelAdjustmentFunction();
                                                                #ifdef CdlEqualFluxinDiffusion
                                                                if(del_dl <= 100.0){
                                                                                                Rdl = ExtraliquidDiffusion/del_dl;
                                                                                                Rp = InterliquidDiffusion/dxCA*2.0*fdP;
                                                                                                CdlCA = (Rp*CdCA + Rdl*CLVFDiffusion)/(Rp + Rdl);
                                                                                }
                                                                                else
                                                                                                CdlCA = CLVFDiffusion;
                                                                #endif

                                                                #ifdef CdlCAConstante
                                                                                CdlCA = CdlCAConst;
                                                                                #endif

                                                                Cdn = CaP.Cd; // Didn't know what to put in this value


                                                                #ifndef ReferenceTest

                                                                #ifdef CdJldExtraTerm
                                                                CdJldVl += -Cdn*(BoundaryConditionHolder[ji].edVd + BoundaryConditionHolder[ji].vcluster);
                                                                #ifndef JldVsDiscrete
                                                                CdJldVs += -Cdn*( BoundaryConditionHolder[ji].vcluster);
                                                                JldVs += -( BoundaryConditionHolder[ji].vcluster);
                                                                #endif
                                                                JldVl += -(BoundaryConditionHolder[ji].edVd + BoundaryConditionHolder[ji].vcluster);

                                                                #endif
                                                                /*
                                                                #ifndef STEREOOneSphereMICRO
                                                                elvlyn[ia*mmax + ja] += dxCA*(-vs + clusters[ClusterP].V[1])/dxVF;
                                                                #endif
                                                                #ifdef STEREOOneSphereMICRO
                                                                elvlyn[ia*mmax + ja] += dxCA*(-vs + clusters[ClusterP].V[1])/dxVF*(M_PI/8.0);
                                                                #endif
                                                                */

                                                                #endif

                                                }

                                                CdConvecTerm = (-BoundaryConditionHolder[ji].edVd)*dt/dyCA*(Cdn-Cds);
                                                #ifdef STEREOOneSphereMICRO
                                                CdConvecTerm = CdConvecTerm*(M_PI/8.0);
                                                #endif
                                                CdCA += CdConvecTerm;
                                //}
                                //else{
                                                //CdConvecTerm = 0.0;
                                //}

                                /*
                                if(clusters[ClusterP].V[1] != 0.0){
                                if(ii == my - 1){
                                        AldOpen += dxCA/dxVF;
                                        AldOpenDown += dxCA/dxVF;
                                }
                                else{
                                if(CaS.state == 0){
                                        AldOpen += dxCA/dxVF;
                                        AldOpenDown += dxCA/dxVF;
                                }

                                }
                                if(ii == 0){
                                        AldOpen -= dxCA/dxVF;
                                        AldOpenUp -= dxCA/dxVF;
                                }
                                else{
                                if(CaN.state == 0){
                                        AldOpen -= dxCA/dxVF;
                                        AldOpenUp -= dxCA/dxVF;
                                }
                                }
                                }
                                */

                                if(BoundaryConditionHolder[ji].on == 1){


                                                if(ii == 0){
                                                        AldOpen -= 1.0;
                                                        AldOpenUp -= 1.0;
                                                        BoundaryConditionHolder[ji].on = 0;
                                                        BoundaryConditionHolder[ji].edVd = 0.0;
                                                        BoundaryConditionHolder[ji].vcluster = 0.0;
                                                        //AldOpenDown -= 1.0;
                                                }
                                                else if(CaN.state == 0.0){
                                                        AldOpen -= 1.0;
                                                        AldOpenUp -= 1.0;
                                                        BoundaryConditionHolder[ji].on = 0;
                                                        BoundaryConditionHolder[ji].edVd = 0.0;
                                                        BoundaryConditionHolder[ji].vcluster = 0.0;
                                                        //AldOpenDown -= 1.0;
                                                }

                                }


                #endif

                #ifndef NOPhaseTransformation
                // Hipótese: Cddt = Ce => TCAdt = Te
                gamaSgamaE = (Teut - TCA - dHCA/soliddensity/Cp)/aTgamaS;
                gamaE = ((Teut - Tf)/ml*(fdCA + afdgamaS*gamaSgamaE) - SUMafdCdvfdCd - afdCdgamaS*gamaSgamaE )/(afdCdgamaE - afdCdgamaS);
                //gamaE = gamaSgamaE;
                if(gamaE > 0.0){
                        gamaS = 0.0;
                        gamaSgamaE = gamaE;
                }
                gamaS = gamaSgamaE - gamaE;
                fECAtdt = fECA - afdgamaS*gamaE;
                fSCAtdt = fSCA - afdgamaS*gamaS;
                fdCAtdt = 1.0 - fSCAtdt - fECAtdt;

                if(fECAtdt > 0.0){
                                ii=ii;
                }

                if(fECAtdt <= 0.0){

                                fECAtdt = 0.0;
                                gamaE = soliddensity/dt*(fECAtdt - fECA);

                                #ifdef ConservativeLocalSolidfraction



                                if(t > 23.277){
                                ii=ii;
                                }

                                a = afdgamaS*aTgamaS;
                                b = aTgamaS*SUMafdvfd + afdgamaS*(dHCA/soliddensity/Cp + TCA - Tf) - ml*afdCdgamaS;
                                c = SUMafdvfd*(dHCA/soliddensity/Cp + TCA - Tf) - ml*SUMafdCdvfdCd;

                                if(a < 0 && c < 0)
                                                ii=ii;

                                if(b*b - 4.0*a*c < 0){
                                                ii=ii;
                                                b = 0.0;
                                                c = 0.0;
                                }

                                gamaS1 = (-b + sqrt(b*b - 4.0*a*c))/2.0/a;
                                gamaS2 = (-b - sqrt(b*b - 4.0*a*c))/2.0/a;

                                TCAtdt = TCA + dHCA/soliddensity/Cp + aTgamaS*gamaS2;
                                if(TCAtdt <= Tf && TCAtdt > Tf + ml*100.0){

                                                gamaS = gamaS2;

                                }

                                TCAtdt = TCA + dHCA/soliddensity/Cp + aTgamaS*gamaS1;
                                if( TCAtdt <= Tf && TCAtdt > Tf + ml*100.0){

                                                gamaS = gamaS1;

                                }
                
                                #endif

                                #ifdef GandinSolidfraction
                                                gamaS = soliddensity/dt*( -dHCA/( soliddensity*Cp*(kpart - 1.0)*(Tliq - Tf)*pow(1.0 - fSCA, kpart - 2.0) + soliddensity*Lf ) );
                                #endif

                                TCAtdt = TCA + dHCA/soliddensity/Cp + aTgamaS*gamaS;
                                fSCAtdt = fSCA - afdgamaS*gamaS;

                                if( fSCAtdt <= 0.0 ){
                                                fSCAtdt = 0.0;
                                                gamaS = soliddensity*(fSCAtdt - 1.0 + SUMafdvfd)/dt;
                                }

                                if( fSCAtdt >= fSmax ){
                                                fSCAtdt = fSmax;
                                                gamaS = soliddensity*(fSCAtdt - 1.0 + SUMafdvfd)/dt;
                                }

                }

                if(fdCAtdt < 1.0 - fSmax){
                                fdCAtdt = 1.0 - fSmax;
                                gamaSgamaE = (fdCAtdt - 1.0 + fSCA + fECA)/afdgamaS;
                                TCAtdt = TCA + dHCA/soliddensity/Cp + aTgamaS*(gamaSgamaE);
                                //gamaE = (fdCAtdt*(Teut - Tf)/ml - SUMafdCdvfdCd - afdCdgamaS*gamaSgamaE)/(afdCdgamaE - afdCdgamaS);
                                gamaE = gamaSgamaE;
                                gamaS = gamaSgamaE - gamaE;
                                fSCAtdt = fSCA - afdgamaS*gamaS;
                                fECAtdt = fECA - afdgamaS*gamaE;
                }




                #ifdef StorageTinCA
                                TCAtdt = TCA + dHCA/soliddensity/Cp + aTgamaS*(gamaS + gamaE);
                                CA[ii*mx + ji].Cddt = (TCAtdt - Tf)/ml;
                #endif

                #ifdef CddtConservative
                CA[ii*mx + ji].Cddt = (SUMafdCdvfdCd + afdCdgamaS*gamaS + afdCdgamaE*gamaE)/(SUMafdvfd + afdgamaS*gamaS + afdgamaS*gamaE);
                #endif

                #ifdef CddtNonConservative
                CA[ii*mx + ji].Cddt = (SUMafdCdvfdCd + afdCdgamaS*gamaS + afdCdgamaE*gamaE)/fdCA ;
                #endif

                #ifdef CddtThermalCAEquilibriumOnlyInFuture
                CA[ii*mx + ji].Cddt = (TCAtdt - Tf)/ml;
                #endif

                #ifdef CddtThermalCAEquilibriumAlways
                CA[ii*mx + ji].Cddt = CA[ii*mx + ji].Cd + (TCAtdt - TCA)/ml;
                #endif

                #ifdef CddtEquilibriumVFt
                CA[ii*mx + ji].Cddt = (VF[ia*mmax + ja].T - Tf)/ml;
                #endif

                #ifdef CASolidConcentration
                if(fSCAtdt != 0.0){
                                CA[ii*mx + ji].Csdt = (CA[ii*mx + ji].fS*CA[ii*mx + ji].Cs - afdCdgamaS*gamaS )/(fSCAtdt);
                }
                else
                                CA[ii*mx + ji].Csdt = CA[ii*mx + ji].Cs;
                #endif

                

                #ifdef CdCAConstante
                        CaP.Cd = CdCAConst;
                        CaP.Cddt = CdCAConst;
                #endif
                #ifdef fSCAConstante
                        CaP.fS = fSCAConst;
                        CaP.fE = 0.0;
                        CaP.dfS = 0.0;
                        CaP.dfE = 0.0;
                #endif


                CA[ii*mx + ji].dfS = fSCAtdt - fSCA;
                CA[ii*mx + ji].dfE = fECAtdt - fECA;

                #endif






                // Para TESTE com um grao!!
                #ifdef NOSolid
                gamaS = 0.0;
                gamaE = 0.0;
                CaP.dfS = 0.0;
                CaP.dfE = 0.0;
                CaP.Cd = Co;
                CaP.Cddt = Co;
                CaP.Csdt = kpart*Co;
                #endif

                #endif // No Scheil

                #ifdef ScheilMicrofS
                /*
                TCAt = getTThermoCouple(t, Ly - (0.5 + ii)*dyCA, ThermoData);
                TCAtdt = getTThermoCouple(t + dt, Ly - (0.5 + ii)*dyCA, ThermoData);

                if(TCAtdt < Tliq)
                CaP.dfS = fSScheil( TCAtdt );
                if(TCA < Tliq)
                        CaP.dfS -= fSScheil(TCAt);

                if(CaP.dfS < 0.0)
                        ii=ii;
                if(CaP.fS + CaP.dfS > fSmax)
                        CaP.dfS = fSmax - CaP.fS;
                */
                #endif

                VF[ia*mmax + ja].dfS += CA[ii*mx + ji].dfS;
                VF[ia*mmax + ja].dfE += CA[ii*mx + ji].dfE;
                clusters[ClusterP].fS += (CaP.dfS + CaP.dfE)*dxCA*dyCA/(clusters[ClusterP].size*dxCA*dyCA);





}
#endif


#ifndef Eutetic
void dfSMicro(){

                if( fabs(CaP.state) >= 2 ){



                                double gamaS;
                                double fSCAtdt;
                                double TCAtdt;
                                #ifdef Eutetic
                                double gamaE;
                                double fECAdt;
                                double fdCAdt;
                                #endif

                                dHCA = dHInterpol();
                                TCA = TInterpol();
                                CdCA = CA[ii*mx + ji].Cd;
                                fSCA = CA[ii*mx + ji].fS;

                                #ifdef eglocal
                                                fGCAdt = (2.0*CaP.L[0]*CaP.L[0]- CaP.Amin)/(CaP.Amax - CaP.Amin);
                                                if(fGCAdt < 0.0)
                                                                fGCAdt = 0.0;
                                                else if(fGCAdt > 1.0)
                                                                fGCAdt = 1.0;
                                                fGCA = fGCAdt - CA[ii*mx + ji].dfG;
                                #endif

                                getCoefficientdfSMicro();




                                #ifdef ConservativeLocalSolidfraction

                                                double a,b,c;
                                                double gamaS1;
                                                double gamaS2;


                                                if(t > 23.277){
                                                ii=ii;
                                                }

                                                a = afdgamaS*aTgamaS;
                                                b = aTgamaS*SUMafdvfd + afdgamaS*(dHCA/soliddensity/Cp + TCA - Tf) - ml*afdCdgamaS;
                                                c = SUMafdvfd*(dHCA/soliddensity/Cp + TCA - Tf) - ml*SUMafdCdvfdCd;

                                                if(a < 0 && c < 0)
                                                                ii=ii;

                                                if(b*b - 4.0*a*c < 0)
                                                                ii=ii;

                                                gamaS1 = (-b + sqrt(b*b - 4.0*a*c))/2.0/a;
                                                gamaS2 = (-b - sqrt(b*b - 4.0*a*c))/2.0/a;

                                                TCAtdt = TCA + dHCA/soliddensity/Cp + aTgamaS*gamaS2;
                                                if(TCAtdt <= Tf && TCAtdt > Tf + ml*100.0){

                                                                gamaS = gamaS2;

                                                }

                                                TCAtdt = TCA + dHCA/soliddensity/Cp + aTgamaS*gamaS1;
                                                if( TCAtdt <= Tf && TCAtdt > Tf + ml*100.0){

                                                                gamaS = gamaS1;

                                                }
                
                                #endif

                                #ifdef GandinSolidfraction
                                                gamaS = soliddensity/dt*( -dHCA/( soliddensity*Cp*(kpart - 1.0)*(Tliq - Tf)*pow(1.0 - fSCA, kpart - 2.0) + soliddensity*Lf ) );
                                #endif

                                TCAtdt = TCA + dHCA/soliddensity/Cp + aTgamaS*gamaS;
                                fSCAtdt = fSCA - afdgamaS*gamaS;

                                if( fSCAtdt <= 0.0 ){
                                                fSCAtdt = 0.0;
                                }
                                if( fSCAtdt >= fSmax ){
                                                fSCAtdt = fSmax;
                                }

                                CA[ii*mx + ji].dfS = fSCAtdt - fSCA;
                                gamaS = soliddensity*(fSCAtdt - 1.0 + SUMafdvfd)/dt;

                                #ifdef CddtConservative
                                CA[ii*mx + ji].Cddt = (SUMafdCdvfdCd + afdCdgamaS*gamaS)/(SUMafdvfd + afdgamaS*gamaS);
                                #endif

                                #ifdef CddtNonConservative
                                CA[ii*mx + ji].Cddt = (SUMafdCdvfdCd + afdCdgamaS*gamaS)/(SUMafdvfd + afdgamaS*gamaS) ;
                                #endif

                                #ifdef CddtThermalCAEquilibriumOnlyInFuture
                                CA[ii*mx + ji].Cddt = (TCAtdt - Tf)/ml;
                                #endif

                                #ifdef CddtThermalCAEquilibriumAlways
                                CA[ii*mx + ji].Cddt = CA[ii*mx + ji].Cd + (TCAtdt - TCA)/ml;
                                #endif

                                #ifdef CddtEquilibriumVFt
                                CA[ii*mx + ji].Cddt = (VF[ia*mmax + ja].T - Tf)/ml;
                                #endif

                                #ifdef CASolidConcentration
                                if(fSCAtdt != 0.0){
                                                CA[ii*mx + ji].Csdt = (CA[ii*mx + ji].fS*CA[ii*mx + ji].Cs - afdCdgamaS*gamaS )/(fSCAtdt);
                                }
                                else
                                                CA[ii*mx + ji].Csdt = CA[ii*mx + ji].Cs;
                                #endif

                                if(CaP.Cddt >= 100.0)
                                                ii=ii;

                                VF[ia*mmax + ja].dfS += CA[ii*mx + ji].dfS;

                }


}



#endif








