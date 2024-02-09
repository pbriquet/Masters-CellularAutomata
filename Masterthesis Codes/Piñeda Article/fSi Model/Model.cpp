

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#pragma hdrstop
#include <condefs.h>

#include "Ahuja.h"
#include "CaMotion.h"
#include "DiffusionLenghts.h"
#include "ExtraFunctions.h"
#include "GlobalDefinitions.h"
#include "GlobalVariables.h"
#include "GrowthKnetics.h"
#include "InterpolationFunctions.h"
#include "NucleationFunctions.h"
#include "OutputFiles.h"
#include "PhaseTransf.h"
#include "SetupDatabase.h"



//#define CpEqScheil

//---------------------------------------------------------------------------
USEUNIT("Ahuja.cpp");
USEUNIT("CaMotion.cpp");
USEUNIT("DiffusionLenghts.cpp");
USEUNIT("ExtraFunctions.cpp");
USEUNIT("GlobalDefinitions.cpp");
USEUNIT("GlobalVariables.cpp");
USEUNIT("GrowthKnetics.cpp");
USEUNIT("InterpolationFunctions.cpp");
USEUNIT("NucleationFunctions.cpp");
USEUNIT("OutputFiles.cpp");
USEUNIT("PhaseTransf.cpp");
USEUNIT("SetupDatabase.cpp");
USEUNIT("SoluteLength.cpp");
//---------------------------------------------------------------------------
#pragma argsused

void SaveInput(){

                FILE * InputFile;
                InputFile = fopen("InputFile","w");



}
void SaveMatrix(){
                FILE * CAFile;
                FILE * VFFile;
                FILE * ClusterFile;

                CAFile = fopen("CAFile","w");
                VFFile = fopen("VFFile","w");
                ClusterFile = fopen("ClusterFile","w");

                fprintf(CAFile, "ii;ji;state;fS;fE;Cd;Cs;L;theta;dC[0];dC[1];dT;Vs;");
                fprintf(VFFile, "ia;ja;fS;fE;fL;CL;T;VL;Vd;Vf;Vs;Se;Ne;");
                fprintf(ClusterFile, "k;thetaclass;cluster;V[0];V[1];d[0];d[1];left;right;up;down;fixed;size;fS;");

                for(RunVFLine){
                for(RunVFCol){
                                fprintf(VFFile, "\n%d;%d;%e;%e;%e;%e;%e;%e;%e;%e;%e;%e;%d;", ia, ja, VfP.fS, VfP.fE, VfP.fL, VfP.CL, VfP.T, VfP.VL, VfP.Vd, VfP.Vf, VfP.Vs, VfP.Se, VfP.Ne);
                for(RunCAinVFLine){
                for(RunCAinVFCol){
                                fprintf(CAFile, "\n%d;%d;%d;%e;%e;%e%e;%e;%e;%e;%e;%e;%e;",ii ,ji , CaP.state , CaP.fS , CaP.fE, CaP.Cd, CaP.Cs, CaP.L[0], CaP.theta, CaP.dC[0], CaP.dC[1], CaP.dT, clusters[ClusterP].V[1]);
                }
                }

                }
                }

                for(RunClusterMatrix){
                                fprintf(ClusterFile, "\n%d;%d;%d;%e;%e;%e;%e;%d;%d;%d;%d;%d;%d;%e;", ii, clusters[ii].thetaclass, clusters[ii].cluster, clusters[ii].V[0], clusters[ii].V[1], clusters[ii].d[0], clusters[ii].d[1], clusters[ii].left, clusters[ii].right, clusters[ii].up, clusters[ii].down, clusters[ii].fixed, clusters[ii].size, clusters[ii].fS);
                }

                fclose(CAFile);
                fclose(VFFile);
                fclose(ClusterFile);

}





//---------------------------------------------------------------------------
int main()
{
                /*
                FILE * AAAHUJAFILE;
                printAhujaFile(AAAHUJAFILE);
                printLevenlspiel();
                */

                FILE *PositionsVF;
                PositionsVF = fopen("PositionVFs.dat", "w");
                for(RunVFLine){
                        fprintf(PositionsVF, "\n%.3f;", Ly - (ia + 0.5)*dyVF);
                }
                fclose(PositionsVF);

                //getTfar();

                double CsVF[mmacroy];
                double CVF[mmacroy][5];
                double tmaxaux;
                double Tout;

                int ClusterMaster;
                bool TCACalculated = 0;

                // Check Compatibility
                if( dxCA != dyCA )
                                printf("\n\n dxCA != dyCA Error.");

                #ifdef NUCLEATIONEqualNinVFs
                                nvmax = NequalinVF;
                #endif

                printf("\nMemory Usage Micro = %d bytes", 2*sizeof(struct MicroMatrix)*mx*my);
                printf("\nMemory Usage Micro = %d bytes", sizeof(struct ClusterMatrix)*mx*my);
                printf("\nMemory Usage Macro = %d bytes", 2*sizeof(struct MacroMatrix)*mmax*mmay);
                printf("\ndxCA = %e", dxCA );
                printf("\ndyCA = %e", dyCA );
                printf("\ndTv = %.2f", dTvnuc );
                printf("\ndTvsigma = %.2f", dTvsigma );
                #ifndef NOMove
                printf("\nMovement Activated");
                #endif


        double tteste = 0.0;        
        InitiateMatrix();
        NucleationStart();
        InitialVisual();

        #ifdef ThermoCouples

        ThermoData = ReadThermoCouple();
        /*
        ja = 0;
        for(RunTime){
                RefreshVFThermoCoupleData(t , ThermoData);

                if(ifProfilesTimeInterval){
                        PrintProfiles();
                }
        }
        return 0;
        */
        #endif
        
        FILE * Jdl;
        FILE * dHFile;

        FILE * CVFall;
        CVFall = fopen("Evo_CVF.dat","w");
        fprintf(CVFall, "t;");
        for(RunVFLine){
                fprintf(CVFall, "VF%d;", ia);
        }
        FILE * fdCdVFall;
        fdCdVFall = fopen("Evo_fdCdVF.dat","w");
        fprintf(fdCdVFall, "t;");
        for(RunVFLine){
                fprintf(fdCdVFall, "VF%d;", ia);
        }
        FILE * flClVFall;
        flClVFall = fopen("Evo_flClVF.dat","w");
        fprintf(flClVFall, "t;");
        for(RunVFLine){
                fprintf(flClVFall, "VF%d;", ia);
        }
        FILE * fsCsVFall;
        fsCsVFall = fopen("Evo_fsCsVF.dat","w");
        fprintf(fsCsVFall, "t;");
        for(RunVFLine){
                fprintf(fsCsVFall, "VF%d;", ia);
        }
        FILE * feCeVFall;
        feCeVFall = fopen("Evo_feCeVF.dat","w");
        fprintf(feCeVFall, "t;");
        for(RunVFLine){
                fprintf(feCeVFall, "VF%d;", ia);
        }

        FILE * TCouple;
        TCouple = fopen("TCouple.dat","w");
        double TCouplePos[3];
        TCouplePos[0] = 0.0263;
        TCouplePos[1] = 0.0603;
        TCouplePos[2] = 0.0933;
        fprintf(TCouple, "t;");
        for(ia = 0; ia < 4; ia++)
                fprintf(TCouple, "TC%d;", ia + 1);


        FILE * DaveFile;
        DaveFile = fopen("Dave.dat","w");
        for(RunClusterMatrix)
                fprintf(DaveFile, "\n%f;",clusters[ii].Dave);

        fclose(DaveFile);

        #ifdef tmaxperV
        tmaxaux = Ly/fabs(clusters[0].V[1]);
        #endif

        dHFile = fopen("dHFile.dat","w");
        fprintf(dHFile, "t;dH;dT;dfS;dHcanculated;");
        //printf("\n\nLf/Cp = %e ; 1.0/pCp = %e", Lf/Cp, 1.0/Cp/soliddensity);
        //system("pause");

        Jdl = fopen("Jdl.dat","w");
        fprintf(Jdl, "t;JdlVF;JDdlVF;CdlJdl");

        //Caverage(CVF);
        AverageCsCdCl(CVF);

        for(RunTime){

                #ifdef ThermoCouples
                RefreshVFThermoCoupleData(t , ThermoData);
                #endif
                
                // Print Results
                ErroCAMotion = 0;
                if(Tboundary( VF[mmacroy - 1].T, Tfarfielddown, Ksolid, Hconvectiondown, dyVF) >= Tliq + 1.0){
                        dt = pow(10.0 , (int) log10(soliddensity*Cp/Ksolid*dyVF*dyVF) - 1.0 );
                }
                else{
                        dt = dtmax;
                }
                #ifdef PRINTTIME
                if(ifTimeTimeInterval){
                        printf("\rt = %4.4f", t);
                }
                #endif

                #ifdef PRINTPROFILES
                if(ifProfilesTimeInterval){
                        PrintProfiles();
                        //Csolidaverage(CsVF);
                        //Caverage(CVF);
                        //PrintCall(CVF);
                        /*
                        fprintf(CVFall, "\n%f;", t);
                        fprintf(fdCdVFall, "\n%f;", t);
                        fprintf(flClVFall, "\n%f;", t);
                        fprintf(fsCsVFall, "\n%f;", t);
                        fprintf(feCeVFall, "\n%f;", t);
                        for(RunVFLine){
                                fprintf(CVFall, "%f;", CVF[ia][0]);
                                fprintf(fdCdVFall, "%f;", CVF[ia][1]);
                                fprintf(flClVFall, "%f;", CVF[ia][2]);
                                fprintf(fsCsVFall, "%f;", CVF[ia][3]);
                                fprintf(feCeVFall, "%f;", CVF[ia][4]);
                        }
                        */
                        //Tout = Tfar();

                        /*
                        fprintf(TCouple, "\n%f;", t);
                        for(ia = 0; ia < 3; ia ++)
                                fprintf(TCouple, "%f;", ThermoCouplePositionT(TCouplePos[ia]) - 273.0);
                        fprintf(TCouple, "%f;", Tout - 273.0);

                        DaveFile = fopen("Dave.dat","w");
                        for(RunClusterMatrix)
                                fprintf(DaveFile, "\n%f;",clusters[ii].Dave);

                        fclose(DaveFile);
                        */


                }
                #endif

                if(t >= 0.222 && t < 0.222 + dt)
                        PrintGrains();

                #ifdef PRINTGRAINS
                if(ifGrainsTimeInterval){

                        PrintGrains();
                        PrintVerticalValues(CVF);

                }
                #endif

                #ifdef PRINTCONTORNO
                if(ifGrainsTimeInterval){
                        PrintContorno();
                }
                #endif

                #ifdef PRINTSURFACECD
                if(ifGrainsTimeInterval){
                                PrintSurfaceCd();
                }
                #endif





                // Calculate dH for VFs


                dHDiffusive();


                #ifdef CpEqScheil


                for(RunVFLine){
                        for(RunVFCol){
                        if( VfP.T <= Tliq && VfP.T > Teut){
                                if(VfP.fS + VfP.fE < fSmax){
                                        VfP.dfS = -VfP.dH/(soliddensity*Cp*(Tliq - Tf)*(kpart - 1.0)*pow(1.0 - VfP.fS , kpart - 2.0) + soliddensity*Lf);

                                }
                                else{
                                        VfP.dfS = 0.0;
                                }



                        }
                        if( VfP.T <= Teut ){
                                VfP.dfS = 0.0;
                                if(VfP.fS + VfP.fE < fSmax){
                                        VfP.dfE = -VfP.dH/soliddensity/Lf;
                                }
                                else
                                        VfP.dfE = 0.0;
                        }

                        dTMacro();
                        }
                }

                RefreshTime();
                #endif

                #ifndef CpEqScheil


                #ifdef fSStateFunction
                for(RunClusterMatrix){
                        #ifndef esiCluster
                        clusters[ii].fS = 0.0;
                        #endif
                        #ifdef esiCluster
                        clusters[ii].fS = clusters[ii].esidt;
                        clusters[ii].countint = 0;
                        #endif
                }
                #endif


                for(RunVFLine){
                        for(RunVFCol){
                                VFP = &VfP;
                                // Calculate dfS and JdlCA for each CA
                                for(RunCAinVFLine){
                                        for(RunCAinVFCol){
                                                CAP = &CaP;
                                                TCACalculated = 0;
                                                //CAP.dT = CAP.dT;

                                                if( CA[ii*mx + ji].dT >= 0.0 || fabs(CA[ii*mx + ji].state) >= 2 ){
                                                        if( ii <= (mmay - 1)*mmiy)
                                                                ii=ii;

                                                }


                                                #ifdef GrowthAlgorithm

                                                if( CA[ii*mx + ji].state >= 2){
                                                        if( fabs(CA[ii*mx + ji].theta) <= M_PI/4.0 + dtheta ){
                                                                //TCA = TInterpol();
                                                                TCA = getTThermoCouple(t, CAPyPosition, ThermoData);
                                                                TCAdt = getTThermoCouple(t + dt, CAPyPosition, ThermoData);
                                                                TCACalculated = 1;
                                                                CAGrowth();
                                                                #ifdef fSStateFunction
                                                                        #ifndef esiCluster
                                                                        #ifdef ScheilMicrofS  // Calcula a fra��o de s�lido local pela equa��o de Scheil
                                                                        CaP.fS = fSScheil(TCA);
                                                                        #endif
                                                                        #ifdef GlobuliticCase
                                                                        CaP.fS = fSmax;
                                                                        #endif
                                                                        #ifdef ConstantMicrofS
                                                                        CaP.fS = ConstantfS;
                                                                        #endif
                                                                        #ifdef LeversRule

                                                                        #endif
                                                                        clusters[ClusterP].fS += CaP.fS;
                                                                        #endif
                                                                #endif

                                                                #ifdef esiCluster
                                                                        if(clusters[ClusterP].fixed == 0.0){
                                                                                clusters[ClusterP].Cdt += (TCA - Tf)/ml;
                                                                                clusters[ClusterP].dCddt += (TCAdt - TCA)/ml/dt;
                                                                                if(CaP.state >= 2){
                                                                                        clusters[ClusterP].ve += DendriticVelocity();
                                                                                        clusters[ClusterP].countint += 1;
                                                                                        /*
                                                                                        if(CaN.state != 0.0)
                                                                                                clusters[ClusterP].Ae += 1.0;
                                                                                        if(CaE.state != 0.0)
                                                                                                clusters[ClusterP].Ae += 1.0;
                                                                                        if(CaS.state != 0.0)
                                                                                                clusters[ClusterP].Ae += 1.0;
                                                                                        if(CaW.state != 0.0)
                                                                                                clusters[ClusterP].Ae += 1.0;
                                                                                        */
                                                                                }

                                                                                //clusters[ClusterP].
                                                                        }
                                                                #endif
                                                        }

                                                }
                                                #endif

                                                #ifdef NucleationAlgorithm
                                                if(VfP.Nuc != 0 && VfP.fL >= dxCA*dxCA/dxVF/dyVF){
                                                        if(TCACalculated == 0)
                                                                TCA = TInterpol();
                                                        Nucleation();
                                                }
                                                #endif
                                        }
                                }
                        }
                }




                for(RunCALine){
                        for(RunCACol){
                                if( fabs(CA[ii*mx + ji].state) >= 2){

                                        ClusterFormation();

                                }
                        }
                }



                #ifndef NoMove

                for(RunClusterMatrix){


                        ClusterMaster = clusters[ii].cluster;

                        while(ClusterMaster != clusters[ClusterMaster].cluster)
                                ClusterMaster = clusters[ClusterMaster].cluster;
                        clusters[ii].cluster = ClusterMaster;

                        if(clusters[ii].cluster == ii ){

                                if(clusters[ii].size < 0)
                                        ii=ii;
                                #ifdef fSStateFunction           // Ajusta o c�lculo de fS da cluster para o tamanho
                                        #ifndef esiCluster
                                        clusters[ii].fS = clusters[ii].fS/clusters[ii].size;
                                        #endif
                                        #ifdef esiCluster
                                        clusters[ii].de = EquivalentD( clusters[ii].size*dxCA*dyCA );
                                        clusters[ii].Ve = M_PI/6.0*pow(clusters[ii].de, 3.0);
                                        clusters[ii].Ae = pow(M_PI, 2.0/3.0)*sqrt(3.0)*pow(clusters[ii].de, 2.0);
                                        if(clusters[ii].countint != 0)
                                                clusters[ii].ve = clusters[ii].ve/clusters[ii].countint;
                                        #endif
                                #endif

                                ClusterVelocity();

                                if(clusters[ii].V[1] < 0)
                                        ii=ii;
                                /*
                                #ifdef MotionStabilityVerification
                                if( modulo(clusters[ii].V) > dxCA/dt)
                                        printf("\n\nAlert! Motion Stability Ruined!");
                                #endif
                                */
                        }
                }

                if(ErroCAMotion){

                                dHDiffusive();
                }





                for(RunClusterMatrix){
                        clusters[ii].Dave += fabs(clusters[clusters[ii].cluster].V[1])*dt;
                        if(clusters[ii].cluster == ii && clusters[ii].fixed == 0){
                                ClusterMovement();

                        }
                }



                CellTranslation();

                for(RunCALine){
                        for(RunCACol){
                                if( fabs(CA[ii*mx + ji].state) >= 2){
                                        ClusterFormation();

                                }
                        }
                }
                /*
                for(RunClusterMatrix){

                                ClusterMaster = clusters[ii].cluster;
                                while(ClusterMaster != clusters[ClusterMaster].cluster)
                                                ClusterMaster = clusters[ClusterMaster].cluster;
                                clusters[ii].cluster = ClusterMaster;
                }
                */
                #endif


                VFTranslation();

                #ifndef NOMove
                //fprintf(efvf, "\n%f;", t);
                //fprintf(elvl, "\n%f;", t);


                // Calculate <vs>s for <vf>f calculation
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
                //Vref = VF[0].Vs;
                //VF[0].Vf = -VF[0].fS*VF[0].Vs;
                /*
                for(RunClusterMatrix){
                        if(clusters[ii].cluster == ii && clusters[ii].fixed == 0){
                                clusters[ii].V[1] = 0.0;
                        }
                }
                */


                //NewCdafterConvection();

                #endif

                for(RunVFLine){
                        for(RunVFCol){
                                VFP = &VfP;
                                stateholder = 0;
                                VfP.Ne = 0;
                                VfP.Nuc = 0;
                                for(RunCAinVFLine){
                                        for(RunCAinVFCol){
                                                // Collect N grains for del_dl calculation
                                                CountGrainsVF();



                                        }
                                }

                               

                        }
                }





                AldOpenUp = 0.0;
                AldOpenDown = 0.0;
                // Calculate dfS for CAs and VF
                for(RunInverseVFLine){
                        ia=ia;
                        for(RunVFCol){
                                VFP = &VfP;

                                CdJdlVF = 0.0;
                                JdlVF = 0.0;
                                JDdlVF = 0.0;
                                dfSVF = 0.0;
                                #ifdef CdJldExtraTerm
                                CdJldVl = 0.0;
                                CdJldVs = 0.0;
                                JldVl = 0.0;
                                JldVs = 0.0;
                                AldOpen = 0.0;

                                #endif

                                if(t >= 0.2905 && t <= 0.2905 + dt)
                                        ia=ia;


                                // Calculate dfS and JdlCA for each CA
                                for(RunInverseCAinVFLine){
                                        for(RunCAinVFCol){
                                                CAP = &CaP;
                                                JDdlCA = 0.0;
                                                JdlCA = 0.0;

                                                dfLMacro();

                                                #ifdef MartoranoLiquidBoundaryCondition
                                                if(CaP.state == 0){
                                                        BoundaryConditionHolder[ji].on = 0;
                                                        BoundaryConditionHolder[ji].edVd = 0.0;
                                                        BoundaryConditionHolder[ji].vcluster = 0.0;
                                                }
                                                #endif

                                                if(CaP.state != 0){
                                                        if(CaP.fS + CaP.fE < fSmax)
                                                                dfSMicro();
                                                        else{
                                                                if(dHInterpol()>0.0)
                                                                        dfSMicro();
                                                        }
                                                }

                                                if(CdlCA != Co)
                                                        ii=ii;
                                                if(CaP.Cd != Co)
                                                        ii=ii;

                                        }
                                }

                                if(elvlyn[1] != 0.0)
                                        ii=ii;

                                // Calculate dfS, CL, T for VF
                                dfSMacro();
                                dTMacro();
                                dCLMacro();

                        }
                }

                /*
                if(ProfilesTimeInterval){
                for(RunVFLine){
                                fprintf(efvf, "%f;", (double) efvfyn[ia]);
                                fprintf(elvl, "%f;", (double) elvlyn[ia]);
                }
                }
                */
        
                RefreshTime();
                #ifdef PRINTCADERESULTS
                if(ifCADEResultsTimeInterval){
                        PrintCADE();

                }
                #endif
                #endif



        }

        //Csolidaverage(CsVF);
        //Caverage(CVF);



        return 0;
}


// Algorithm sequence:
        // FOR(t)
        // * return dt to dtmax
        // * print profiles and grains
        // 2) Translation of Grains
        //
        //      # Scanning MicroA:
        //      2.1) Cluster Formations (Read MicroA. If a active CA has a boundary in common with the domain or another active CA, fix the cluster in the first case, and in the second case, gather both cluster into one cluster)
        //
        //      # Scanning ClusterData: (only Grainnumber=Clusternumber, and until Grainnumber = NTotalGrains )
        //      2.2) Calculate Velocities (LaGrangian way.. because just there isn't a way to commmunicate Navier-Stokes results with CA yet)
        //              * get dt for stability
        //      2.3) Movement. Add to vector displacement.
        //
        //      (To Decide)
        //      # Scanning Macro
        //              Get values of dHVF
        //
        //      # Scanning MicroA again:
        //      2.4) Nucleation and CA Growth.
        //
        //      # Scanning MicroA again:
        //      2.5) Copy MicroA to MicroB after translation.
        //
        // 3) Phase Transformation
        //
        //      (To Decide)
        //      # Scanning Macro
        //      3.1) Get values for dHVF(t). (This can be made before without issues.. maybe it is better to get it before copying MicroA into MicroB)
        //
        //      # Scanning Macro
        //      ## Scanning just the MicroB CAs in Macro VF
        //      3.2) Calculate CA local solidification gamaS. Account for JDdl in VF. Already sum in local solid fraction the phase transformation.
        //      3.3) Activate CAs, count as phase transformation in VF.
        //
        // 4) Thermal Balance
        //      4.1) Use sum of gamaS to calculate later T of VF
        //      4.2) Use JDdl and Jgamadl to calculate new liquid concentration.
        //
        // 5) Refresh Matrix
        //      5.1) Copy MicroB to MicroA
        //      5.2) Refresh Macro.
        //      5.3) Refresh time. Loop.
        // END FOR(t)

