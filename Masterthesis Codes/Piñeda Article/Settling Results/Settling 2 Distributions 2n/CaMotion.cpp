//---------------------------------------------------------------------------
#pragma hdrstop

#include <math.h>
#include <stdlib.h>

#include "CaMotion.h"
//---------------------------------------------------------------------------
#pragma package(smart_init)


#include "GlobalDefinitions.h"
#include "SetupDatabase.h"
#include "GlobalVariables.h"

#include "Ahuja.h"


double SoluteDiffusionConvecLength( double El, double RelativeVl, double Ne ){
        double result;

        double Rcell,Re,Eg;

        double Cepsilon, Pecletepsilon;

                if(Ne == 0.0)
                                return 1.0;

        Eg = 1 - El;

        if (Eg > 0.97 ){
                Eg = 0.97;
                El = 1.0 - Eg;
        }
        else if (Eg < 1E-12){
                Eg = 1E-12;
                El = 1.0 - Eg;
        }

        Rcell = pow(3.0/(4.0*M_PI*Ne),1.0/3.0);
        Re = Rcell*pow(Eg,(1.0/3.0));

        Cepsilon = (2.0 + 4.0/3.0*Eg)/(2.0 - 3.0*pow(Eg, 1.0/3.0) + 3.0*pow(Eg, 5.0/3.0) - 2.0*pow(Eg, 2.0));
        Pecletepsilon = El*RelativeVl*2.0*Re/DiffLiq;

        result = 2.0*Re/(2.0 + 0.865*pow( Cepsilon/El ,1.0/3.0)*pow( Pecletepsilon,1.0/3.0));

        return result;
}

        // 2.8) Funções para Cálculo de Velocidade

double fKvpartition(double el, double beta, double betad){
                return (1.0 - el)*pow(beta/betad, 2.0);
}


double fbetad(double esi, double phie, double Senv, double Ssol){
                double xresp = 3.0*sqrt(5.0)*Ssol/Senv/phie/pow(1.0 - esi, 3.0/2.0);
                return  xresp;
}

double fCpphie(double phie, double el){
                if(el >= 0.7){
                                return 1.26*log10(phie/0.163);
                }
                else{
                                return phie*phie;
                }
}

double fbetal(double betad, double el, double phie){
                double eta = pow(1.0 - el, 1.0/3.0);
                double Csepsilon = 1.0/fCpphie(phie, el)*(2.0*pow(betad,2.0)*(1.0 - tanh(betad)/betad ) )/(2.0*pow(betad,2.0) + 3.0*(1.0 - tanh(betad)/betad));
                double Cmepsilon = (2.0 + 4.0/3.0*pow(eta, 5.0) )/(2.0 - 3.0*eta + 3.0*pow(eta,5.0) - 2.0*pow(eta,6.0));
                return pow( 9.0/2.0*(1.0 - el)*Cmepsilon*Csepsilon, 1.0/2.0);
}

double fbeta(double betad, double betal, double el){
                double nbeta = 0.176*log10(betad) + 0.275;
                return betad/pow( pow(1.0 - el, nbeta) + pow(betad/betal, 2.0*nbeta), 1.0/2.0/nbeta);
}

void CalculateAverageVs(){
        if( fabs(CaP.state) >= 2){
                if(VfP.fS != 0.0)
                        ii=ii;
                //VfP.Vs += VfP.fL!=1.0 ? (CaP.fS != 0.0 ? clusters[ClusterP].V[1]/(1.0 - VfP.fL)/mmix/mmiy : 0.0) : 0.0;
                if(VfP.fS !=0.0){
                        if(CaP.fS + CaP.fE != 0.0)
                                VfP.Vs += (CaP.fS + CaP.fE)*clusters[ClusterP].V[1]/(VfP.fS + VfP.fE)/mmiy/mmix;
                }

        }

}

void CalculateAverageVfVdVl(){
        double Kv, betad, betal, beta, Senv, Ssol, el, esi, phie,ef,ed;

        if(VfP.Vs != 0.0){
                phie = 0.846;

                //el = VfP.fL;
                el = VfP.fLCA;
                ef = 1.0 - VfP.fS - VfP.fE;
                ed = ef - el;
                if(el == 1.0){
                        VfP.Vf = 0.0;
                        VfP.Vd = 0.0;
                        VfP.VL = 0.0;
                }
                else{
                        esi = VfP.fS/(1.0 - el);

                        Senv = VfP.Se*(4.0/M_PI);
                        //Senv = 625.0;
                        Ssol = 2.0/lambda2;
                        //Ssol = 0.1*Senv;

                        betad = fbetad( esi, phie , Senv, Ssol );
                        betal = fbetal(betad,el,phie);
                        beta = fbeta(betad,betal,el);
                        #ifdef efVfPequalseSvSP
                        VfP.Vf = -(VfP.fS + VfP.fE)*VfP.Vs/ef;
                        #endif
                        #ifdef efVfPequalseSvSN
                        if(ia != 0)
                        VfP.Vf = -(VfN.fS + VfN.fE)*VfN.Vs/ef;
                        else
                        VfP.Vf = 0.0;
                        #endif

                        if(1.0 - VfP.fS - VfP.fL > fdminVF){
                                Kv = fKvpartition(el,beta,betad);
                                if(Kv >= 0.001)
                                        ii=ii;
                                VfP.Vd = VfP.Vs + Kv*ef/ed*(VfP.Vf - VfP.Vs);
                                VfP.VL = VfP.Vs + (1.0 - Kv)*ef/el*(VfP.Vf - VfP.Vs);
                        }
                        else{
                                VfP.Vd = VfP.Vs;
                                VfP.VL = VfP.Vf;
                        }
                }
        }
        else{
        VfP.Vf = 0.0;
        VfP.Vd = 0.0;
        VfP.VL = 0.0;
        }

}





double EquivalentD( double area ){

        double  result;

        #ifdef SphereCA
        result = 2*pow( area/M_PI , 0.5);
        #endif
        #ifdef OctahedronCA
        double dc = sqrt(area);

        result = sqrt( 2.0 )*dc/pow( M_PI, 1.0/3.0 );
        #endif
        return result;
}


double DimensionlessDiameter(double dsph, double clusterdensity){
        double result;
        result = dsph*pow( gravity*liquiddensity*fabs((clusterdensity - liquiddensity))/pow(viscosity, 2.0)     , 1.0/3.0);
        return result;
}

double TerminalVelocity( double DimensionlessV, double clusterdensity){
        double result;
        result = DimensionlessV*pow( pow(liquiddensity, 2.0)/( gravity*viscosity*fabs((clusterdensity - liquiddensity)) )    , -1.0/3.0);
        return result;
}

double DimensionlessVelocity( double DimensionlessD){
        double result;
        result = pow( pow( 18.0/pow(DimensionlessD ,2.0) , K2) + pow( 3*K1/4.0/pow( DimensionlessD  , 0.5)  , K2)   , -1.0/K1 );
        return result;
}


double CalculateVelocity( double area , double clusterdensity ){
        double result =0.0;
        if(area != 0.0 && clusterdensity != liquiddensity){
                result = TerminalVelocity( DimensionlessVelocity( DimensionlessDiameter( EquivalentD( area ), clusterdensity ) ), clusterdensity );

                if(clusterdensity < liquiddensity)
                        result = -result;
        }
        else
                result = 0.0;

        if(fabs(result) > dyCA/dt){
                //printf("\nMotion Criteria in Calculate Velocity");
                //system("pause");
                dt = dyCA/fabs(result)*0.5;
                ErroCAMotion = 1;
        
        }


        //result = 0.0;
        return result;
}


void ClusterFormation(){

        int ClusterPoint = ClusterP;
        int ClusterMaster;
        int ClusterEastNeighbor;
        int ClusterSouthNeighbor;
        //int ClusterNorthNeighbor;
        //int ClusterWestNeighbor;

        int iaux;
        iaux = ClusterPoint;
        while(clusters[iaux].cluster != iaux){
                clusters[iaux].cluster = min(clusters[iaux].cluster,iaux);
                iaux = min(clusters[iaux].cluster,iaux);
        }
        ClusterP = iaux;
        ClusterPoint = iaux;

        ClusterMaster = ClusterPoint;
        ClusterEastNeighbor = ClusterPoint;
        ClusterSouthNeighbor = ClusterPoint;
        //ClusterNorthNeighbor = ClusterPoint;
        //ClusterWestNeighbor = ClusterPoint;

        if( ji != mx - 1 && fabs(CaE.state) >= 2){
                iaux = ClusterE;
                while(clusters[iaux].cluster != iaux){
                        clusters[iaux].cluster = min(clusters[iaux].cluster,iaux);
                        iaux = min(clusters[iaux].cluster,iaux);
                }
                ClusterE = iaux;
                ClusterEastNeighbor = ClusterE;

                }
        else{
                ClusterEastNeighbor = ClusterPoint;
        }

        if(ii != my - 1 && fabs(CaS.state) >= 2){
                iaux = ClusterS;
                while(clusters[iaux].cluster != iaux){
                        clusters[iaux].cluster = min(clusters[iaux].cluster,iaux);
                        iaux = min(clusters[iaux].cluster,iaux);
                }
                ClusterS = iaux;
                ClusterSouthNeighbor = ClusterS;
        }
        else{
                ClusterSouthNeighbor = ClusterPoint;
        }

        /*

        if(ii != 0 && fabs(CaN.state) >= 2){
                iaux = ClusterN;
                while(clusters[iaux].cluster != iaux){
                        clusters[iaux].cluster = min(clusters[iaux].cluster,iaux);
                        iaux = min(clusters[iaux].cluster,iaux);
                }
                ClusterN = iaux;
                ClusterNorthNeighbor = ClusterN;
        }
        else{
                ClusterNorthNeighbor = ClusterPoint;
        }

        if(ji != 0 && fabs(CaW.state) >= 2){
                iaux = ClusterW;
                while(clusters[iaux].cluster != iaux){
                        clusters[iaux].cluster = min(clusters[iaux].cluster,iaux);
                        iaux = min(clusters[iaux].cluster,iaux);
                }
                ClusterW = iaux;
                ClusterWestNeighbor = ClusterW;
        }
        else{
                ClusterWestNeighbor = ClusterPoint;
        }
        */


        // Fixed Verification

        if(clusters[ClusterPoint].fixed == 0 && (FixClusterConditions) ){
                clusters[ ClusterPoint ].fixed = 1;
                clusters[ ClusterPoint ].V[0] = 0.0;
                clusters[ ClusterPoint ].V[1] = 0.0;
        }




        // Right Verification

        if( ClusterEastNeighbor != ClusterPoint){

                if(ClusterEastNeighbor < ClusterPoint )
                        ClusterMaster = ClusterEastNeighbor;
                else
                        ClusterMaster = ClusterPoint ;

                if( clusters[ClusterPoint ].size + clusters[ClusterEastNeighbor].size > mx*my)
                        ii=ii;

                clusters[ClusterMaster].fS = (clusters[ClusterPoint ].size*clusters[ClusterPoint ].fS + clusters[ClusterEastNeighbor].size*clusters[ClusterEastNeighbor].fS)/(clusters[ClusterPoint].size + clusters[ClusterEastNeighbor].size);
                clusters[ClusterMaster].size = clusters[ClusterPoint ].size + clusters[ClusterEastNeighbor].size;

                if( clusters[ClusterPoint].fixed == 1 || clusters[ClusterEastNeighbor].fixed == 1){
                        clusters[ClusterMaster].fixed = 1;
                        clusters[ClusterMaster].V[1] = 0.0;
                }

                clusters[ClusterE].cluster = ClusterMaster;
                clusters[ClusterP].cluster = ClusterMaster;
                ClusterE = ClusterMaster;
                ClusterP = ClusterMaster;


        }

        // Down Verification

        if( ClusterSouthNeighbor != ClusterPoint ){

                if(ClusterSouthNeighbor < ClusterPoint )
                        ClusterMaster = ClusterSouthNeighbor;
                else
                        ClusterMaster = ClusterPoint ;

                if( clusters[ClusterPoint ].size + clusters[ClusterSouthNeighbor].size > mx*my)
                        ii=ii;

                clusters[ClusterMaster].fS = (clusters[ClusterPoint ].size*clusters[ClusterPoint ].fS + clusters[ClusterSouthNeighbor].size*clusters[ClusterSouthNeighbor].fS)/(clusters[ClusterPoint ].size + clusters[ClusterSouthNeighbor].size);
                clusters[ClusterMaster].size = clusters[ClusterPoint ].size + clusters[ClusterSouthNeighbor].size;

                if( clusters[ClusterPoint ].fixed == 1 || clusters[ClusterSouthNeighbor].fixed == 1){
                        clusters[ClusterMaster].fixed = 1;
                        clusters[ClusterMaster].V[1] = 0.0;
                }


                clusters[ClusterS].cluster = ClusterMaster;
                clusters[ClusterP].cluster = ClusterMaster;
                ClusterS = ClusterMaster;
                ClusterP = ClusterMaster;

        }



        /*
        if( ClusterNorthNeighbor != ClusterPoint){

                if(ClusterNorthNeighbor < ClusterPoint)
                        ClusterMaster = ClusterNorthNeighbor;
                else
                        ClusterMaster = ClusterPoint ;

                clusters[ClusterMaster].fS = (clusters[ClusterPoint ].size*clusters[ClusterPoint ].fS + clusters[ClusterNorthNeighbor].size*clusters[ClusterNorthNeighbor].fS)/(clusters[ClusterPoint ].size + clusters[ClusterNorthNeighbor].size);
                clusters[ClusterMaster].size = clusters[ClusterPoint ].size + clusters[ClusterNorthNeighbor].size;

                if( clusters[ClusterPoint].fixed == 1 || clusters[ClusterNorthNeighbor].fixed == 1)
                        clusters[ClusterMaster].fixed = 1;

                ClusterN = ClusterMaster;
                ClusterP = ClusterMaster;

        }

        if( ClusterWestNeighbor != ClusterPoint){

                if(ClusterWestNeighbor < ClusterPoint)
                        ClusterMaster = ClusterWestNeighbor;
                else
                        ClusterMaster = ClusterPoint ;

                clusters[ClusterMaster].fS = (clusters[ClusterPoint ].size*clusters[ClusterPoint ].fS + clusters[ClusterWestNeighbor].size*clusters[ClusterWestNeighbor].fS)/(clusters[ClusterPoint ].size + clusters[ClusterWestNeighbor].size);
                clusters[ClusterMaster].size = clusters[ClusterPoint ].size + clusters[ClusterWestNeighbor].size;

                if( clusters[ClusterPoint].fixed == 1 || clusters[ClusterWestNeighbor].fixed == 1)
                        clusters[ClusterMaster].fixed = 1;

                ClusterW = ClusterMaster;
                ClusterP = ClusterMaster;

        }
        */
        if(clusters[ClusterMaster].size > mx*my)
                        ii=ii;


        /*
        if(clusters[ClusterMaster].cluster!=ClusterMaster){
                int ClusterMasterI = (clusters[ClusterP].cluster);
                int ClusterMasterII = clusters[ClusterMasterI].cluster;
                int ClusterSlave = ClusterP;

                ClusterP=ClusterMasterII;
        }
        */


}



void ClusterVelocity(){

        #ifdef MotionAlgorithm
        #ifdef LevenspielVelocity
        if(clusters[ii].fixed == 0)
        clusters[ii].V[1] = -CalculateVelocity( clusters[ii].size*dxCA*dyCA , clusters[ii].fS*soliddensity + (1.0 - clusters[ii].fS)*liquiddensity );
        #endif

        #ifdef ConstantTerminalVelocity
        if(clusters[ii].fixed == 0)
        clusters[ii].V[1] = -0.01;
        #endif
        #ifdef AhujaVelocity
        double LevenspielComparison = 0.0;
        if(clusters[ii].fixed == 0){
        clusters[ii].V[1] = -CalculateVClusterAhuja();

        LevenspielComparison = -CalculateVelocity( clusters[ii].size*dxCA*dyCA , clusters[ii].fS*soliddensity + (1.0 - clusters[ii].fS)*liquiddensity );
        if(LevenspielComparison != clusters[ii].V[1])
                ii=ii;
        }
        #endif
        else
                clusters[ii].V[1] = 0.0;
        #endif
        #ifndef MotionAlgorithm
                clusters[ii].V[1] = 0.0;
        #endif

        if(fabs(clusters[ii].V[1])>dyCA/dt){
                printf("\n\nWarning. CA Motion Criteria");
        }

}


void ClusterMovement(){

        clusters[ii].d[0] += clusters[ii].V[0]*dt;
        clusters[ii].d[1] += clusters[ii].V[1]*dt;

        if( clusters[ii].d[0] > dxCA/2.0){
                clusters[ii].right = 1;
                clusters[ii].d[0] -= dxCA;
        }

        if( clusters[ii].d[0] < -dxCA/2.0){
                clusters[ii].left = 1;
                clusters[ii].d[0] += dxCA;
        }

        if( clusters[ii].d[1] > dyCA/2.0){
                clusters[ii].up = 1;
                clusters[ii].d[1] -= dyCA;
        }

        if( clusters[ii].d[1] < -dyCA/2.0){
                clusters[ii].down = 1;
                clusters[ii].d[1] += dyCA;
        }


}


void CellTranslation(){

                for(RunVFLine){
                for(RunVFCol){
                for(RunCAinVFLine){
                                for(RunCAinVFCol){
                                                CaPaux.state = CaP.state;
                                                CaPaux.fS = CaP.fS;
                                                CaPaux.dfS = CaP.dfS;
                                                CaPaux.theta = CaP.theta;
                                                CaPaux.L[0] = CaP.L[0];
                                                CaPaux.L[1] = CaP.L[1];
                                                CaPaux.L[2] = CaP.L[2];
                                                CaPaux.L[3] = CaP.L[3];
                                                CaPaux.L[4] = CaP.L[4];
                                                CaPaux.dC[0] = CaP.dC[0];
                                                CaPaux.dC[1] = CaP.dC[1];
                                                CaPaux.Cd = CaP.Cd;
                                                CaPaux.Cddt = CaP.Cddt;
                                                CaPaux.Cs = CaP.Cs;
                                                CaPaux.Csdt = CaP.Csdt;

                                                #ifdef KeepClinCA
                                                CaPaux.Cl = CaP.Cl;
                                                #endif

                                                #ifdef eglocal
                                                CaPaux.Amin = CaP.Amin;
                                                CaPaux.Amax = CaP.Amax;
                                                #endif

                                                #ifdef Eutetic
                                                CaPaux.fE = CaP.fE;
                                                CaPaux.dfE = CaP.dfE;
                                                #endif

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
                                                CaP.Cd = VfP.CL;
                                                CaP.Cs = 0.0;
                                                CaP.Csdt = 0.0;

                                                #ifdef eglocal
                                                CaP.Amin = 0.0;
                                                CaP.Amax = 0.0;
                                                #endif

                                                #ifdef Eutetic
                                                CaP.fE = 0.0;
                                                CaP.dfE = 0.0;
                                                #endif

                                }
                }
                }
                }

                for(RunVFLine){
                for(RunVFCol){
                for(RunCAinVFLine){
                                for(RunCAinVFCol){

                                if( fabs(CaPaux.state) >= 2.0){

                                                if(FuturePosition != ii*mx + ji)
                                                                ii=ii;

                                                CA[FuturePosition].state = CAaux[ii*mx + ji].state;
                                                CA[FuturePosition].fS = CAaux[ii*mx + ji].fS;
                                                CA[FuturePosition].dfS = CAaux[ii*mx + ji].dfS;
                                                CA[FuturePosition].theta = CAaux[ii*mx + ji].theta;
                                                CA[FuturePosition].L[0] = CAaux[ii*mx + ji].L[0];
                                                CA[FuturePosition].L[1] = CAaux[ii*mx + ji].L[1];
                                                CA[FuturePosition].L[2] = CAaux[ii*mx + ji].L[2];
                                                CA[FuturePosition].L[3] = CAaux[ii*mx + ji].L[3];
                                                CA[FuturePosition].dC[0] = CAaux[ii*mx + ji].dC[0];
                                                CA[FuturePosition].dC[1] = CAaux[ii*mx + ji].dC[1];
                                                CA[FuturePosition].Cd = CAaux[ii*mx + ji].Cd;
                                                CA[FuturePosition].Cddt = CAaux[ii*mx + ji].Cddt;
                                                CA[FuturePosition].Cs = CAaux[ii*mx + ji].Cs;
                                                CA[FuturePosition].Csdt = CAaux[ii*mx + ji].Csdt;

                                                #ifdef eglocal
                                                CA[FuturePosition].Amin = CAaux[ii*mx + ji].Amin;
                                                CA[FuturePosition].Amax = CAaux[ii*mx + ji].Amax;
                                                #endif

                                                #ifdef Eutetic
                                                CA[FuturePosition].fE = CAaux[ii*mx + ji].fE;
                                                CA[FuturePosition].dfE = CAaux[ii*mx + ji].dfE;
                                                #endif

                                                
                                }
                                

                                }
                }
                }
                }

}

void VFTranslation(){
                #ifdef MotionAlgorithm
                for(RunCACol){
                                BoundaryConditionHolder[ji].edVd = 0.0;
                                BoundaryConditionHolder[ji].on = 0;
                }
                #endif

                for(RunInverseVFLine){
                                for(RunVFCol){

                                                VFaux[ia*mmax + ja].fS = VfP.fS;
                                                VFaux[ia*mmax + ja].fE = VfP.fE;
                                                VFaux[ia*mmax + ja].fL = VfP.fL;
                                                VfP.fS = 0.0;
                                                VfP.fE = 0.0;

                                                VfP.fLCA = 1.0;

                                                efvfyn[ia*mmax + ja] = 0.0;
                                                elvlyn[ia*mmax + ja] = 0.0;

                                                VfP.Se = 0.0;



                                                for(RunInverseCAinVFLine){
                                                                for(RunCAinVFCol){


                                                                                if( fabs(CaP.state) >= 2.0){



                                                                                                VfP.fS += CaP.fS/mmix/mmiy;
                                                                                                VfP.fE += CaP.fE/mmix/mmiy;
                                                                                                #ifndef eglocal
                                                                                                VfP.fLCA -= 1.0/mmix/mmiy;
                                                                                                #endif
                                                                                                #ifdef eglocal
                                                                                                if( CaP.Amax != 0.0 )
                                                                                                VfP.fLCA -= (2.0*CaP.L[0]*CaP.L[0] - CaP.Amin)/(CaP.Amax - CaP.Amin)/mmix/mmiy;
                                                                                                else
                                                                                                VfP.fLCA -= 0.0;
                                                                                                #endif
                                                                                                
                                                                                                // Turn off Up/Down/Left/Right in Clusters

                                                                                               clusters[ClusterP].right = 0;
                                                                                               clusters[ClusterP].up = 0;
                                                                                               clusters[ClusterP].left = 0;
                                                                                               clusters[ClusterP].down = 0;



                                                                                               if(ii != 0 && CaN.state == 0){
                                                                                                        VfP.Se+=dxCA/dxVF/dyVF;

                                                                                                        
                                                                                                }
                                                                                                if(ii != my - 1 && CaS.state == 0){
                                                                                                        VfP.Se+=dxCA/dxVF/dyVF;


                                                                                                }
                                                                                                if(ji != 0 && CaW.state == 0){
                                                                                                        VfP.Se+=dyCA/dxVF/dyVF;
                                                                                                        }
                                                                                                if(ji != mx - 1 && CaE.state == 0){
                                                                                                        VfP.Se+=dyCA/dxVF/dyVF;
                                                                                                }




                                                                                }
                                                                                else{
                                                                                                #ifdef MotionAlgorithm
                                                                                               BoundaryConditionHolder[ji].on = 0;
                                                                                               #endif
                                                                                }
                                                                }
                                                }

                                                // efvfyn Calculation
                                                if(ia == mmay - 1)
                                                                efvfyn[ia] = (VfP.fS - VFaux[ia*mmax + ja].fS)*dyVF/dt;
                                                else
                                                                efvfyn[ia] = (VfP.fS - VFaux[ia*mmax + ja].fS)*dyVF/dt + efvfyn[ia+1];
                                                if(ia == 0)
                                                        efvfyn[ia] = 0.0;

                                                /*
                                                if(ia == mmay - 1)
                                                                elvlyn[ia] = elvlyn[ia];
                                                else
                                                                elvlyn[ia] += elvlyn[ia+1];
                                                if(ia == 0)
                                                        elvlyn[ia] = 0.0;
                                                */

                                                if( VfP.fLCA < fLVFmin )
                                                                VfP.fLCA = fLVFmin;

                                                //VfP.fL = VfP.fLCA;



                                }
                }

}





void NewCdafterConvection(){

                double edPVdyn, edSVdys;
                double edP;
                double CdS;

                for(RunCALine){
                                for(RunCACol){
                                                CaPaux.Cddt = CaP.Cd;

                                                if(fabs(CaP.state) >= 2){


                                                                if(ii!=0){
                                                                                if((CaNaux.fS + CaNaux.fE) - (CaPaux.fS + CaPaux.fE) != 0.0){
                                                                                                edPVdyn = ( (CaP.fS - CaPaux.fS) + (CaP.fE - CaPaux.fE) )*(CaNaux.fS + CaNaux.fE)/( (CaNaux.fS + CaNaux.fE) - (CaPaux.fS + CaPaux.fE) )*dyCA/dt;
                                                                                }
                                                                                else{
                                                                                                if((CaP.fS - CaPaux.fS) + (CaP.fE - CaPaux.fE) == 0.0){
                                                                                                                edPVdyn = (CaNaux.fS + CaNaux.fE)*dyCA/dt;
                                                                                                }
                                                                                                else{
                                                                                                                edPVdyn = 0.0;
                                                                                                }
                                                                                }

                                                                }
                                                                else{
                                                                                edPVdyn = 0.0;
                                                                }
                                                                if(ii!=my -1){
                                                                                if( (CaP.fS + CaP.fE) - (CaS.fS + CaS.fE) != 0.0){
                                                                                                edSVdys = ( (CaP.fS - CaPaux.fS) + (CaP.fE - CaPaux.fE) )*(CaPaux.fS + CaPaux.fE)/( (CaP.fS + CaP.fE) - (CaS.fS + CaS.fE) )*dyCA/dt;
                                                                                }
                                                                                else{
                                                                                                if( (CaP.fS - CaPaux.fS) + (CaP.fE - CaPaux.fE) == 0.0 ){
                                                                                                                edSVdys = (CaPaux.fS + CaPaux.fE)*dyCA/dt;
                                                                                                }
                                                                                                else{
                                                                                                                edSVdys = 0.0;
                                                                                                }
                                                                                }
                                                                                if(fabs(CaSaux.state) >= 2){
                                                                                                CdS = CaSaux.Cd;
                                                                                }
                                                                                else{
                                                                                                CdS = CaPaux.Cd;
                                                                                }
                                                                }
                                                                else{
                                                                                edSVdys = 0.0;
                                                                                CdS = 0.0;
                                                                }

                                                                if( edPVdyn != 0.0 || edSVdys != 0.0){
                                                                                ii=ii;
                                                                }

                                                                edP = 1.0 - CaPaux.fS - CaPaux.fE;
                                                                /*
                                                                edPVdyn = edP*0.01;

                                                                if(ii != my)
                                                                                edSVdys = edP*0.01;
                                                                else
                                                                                edSVdys = 0.0;
                                                                */
                                                                CaP.Cd = (edP*CaPaux.Cd - (edPVdyn*CaPaux.Cd - edSVdys*CdS)*dt/dyCA)/(1.0 - CaP.fS - CaP.fE);

                                                                ii = ii;

                                                                if(ii != 0.0 && CaPaux.fS + CaNaux.fS >= 1.0){
                                                                                ii=ii;
                                                                }
                                                                if(CaP.Cd < Co)
                                                                                ii=ii;


                                                }
                                }
                }

}



void printLevenlspiel(){
                FILE * LevenspielFile;
                LevenspielFile = fopen("Levenspiel.dat","w");

                fprintf(LevenspielFile, "dstar;ustar");

                double DD;

                for( DD = 1.0; DD < 1.0e4 ; DD += DD){
                                fprintf(LevenspielFile, "\n%.12f;%.12f;", log10(DD), log10(DimensionlessVelocity( DD)));
                }

                fclose(LevenspielFile);
}

