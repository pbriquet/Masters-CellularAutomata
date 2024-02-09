//---------------------------------------------------------------------------
#pragma hdrstop

#include <stdio.h>
#include <math.h>
#include "GlobalDefinitions.h"
#include "SetupDatabase.h"
#include "GlobalVariables.h"

#include "ExtraFunctions.h"
#include "InterpolationFunctions.h"

#include "OutputFiles.h"
//---------------------------------------------------------------------------
#pragma package(smart_init)

// OutputFiles








#ifdef THERMALPROFILE
FILE *TMacro;
#endif

#ifdef SOLIDFRACTIONPROFILE
FILE *fSMacro;
#endif

#ifdef Eutetic
#ifdef EUTETICFRACTIONPROFILE
FILE *fEMacro;
#endif
#ifdef SUMSOLIDEUTETICFRACTIONPROFILE
FILE *fSfEMacro;
#endif
#endif

#ifdef THERMALHOPPROFILE
FILE *TMacroHop;
#endif

#ifdef SOLIDFRACTIONHOPPROFILE
FILE *fSMacroHop;
#endif

#ifdef PRINTGRAINS
FILE *gFile;
#endif

#ifdef PRINTCONTORNO
FILE *cFile;
#endif


#ifdef PRINTSOLUTEPROFILE
FILE *flClFile;
#endif

#ifdef PRINTLIQUIDFRACTIONPROFILE
FILE *fLMacro;
FILE * fLCAMacro;
#endif

#ifdef PRINTAverageVelocities
FILE *VLMacro;
FILE *VdMacro;
FILE *VsMacro;
FILE *elVLMacro;
#endif

#ifdef PRINTSeProfile
FILE *Evolution_Se;
FILE *Evolution_Ne;
FILE *Evolution_Nuc;
#endif

#ifdef PRINTSURFACECD
FILE *SurfaceCd;
FILE *SurfaceVd;
FILE *SurfacefS;
#endif

FILE *efvf;
FILE *elvl;

#ifdef PRINTVERTICALVALUES
FILE *VertCL;
FILE *VertVL;
FILE *VertVs;
FILE *VertVd;
FILE *VertfS;
FILE *VertflVl;
FILE *VertfsVs;
FILE *VertfdVd;
FILE *VertSe;
FILE *VertT;
FILE *VertCs;
FILE *VertC;
#endif

#ifdef PRINTCALL
FILE *Evolution_Call;
#endif


void Csolidaverage(double CsVF[mmacroy]){
                /*
                double esCs,eeCe;
                double Ceut = (Teut - Tf)/ml;
                FILE * CsFile;
                CsFile = fopen("CsFile.dat","w");
                fprintf(CsFile, "y;Cs;");

                for(RunVFLine){
                        for(RunVFCol){
                                for(RunCAinVFLine){
                                        for(RunCAinVFCol){

                                                if(CaP.state != 0){
                                                                esCs += CaP.fS*CaP.Cs;
                                                                eeCe += CaP.fE*Ceut;
                                                }

                                        }
                                }
                                CsVF[ia] = (esCs + eeCe)/mmiy/mmix;
                                esCs = 0.0;
                                eeCe = 0.0;
                                fprintf(CsFile, "\n%f;%f", Ly - dyVF*(ia + 0.5), CsVF[ia] );
                        }

                }
                fclose(CsFile);

                */
}

void Caverage(double CVF[mmacroy]){
                /*
                double esCs,eeCe,edCd,elCl;
                double Ceut = (Teut - Tf)/ml;

                FILE * CFile;
                CFile = fopen("CFile.dat","w");
                fprintf(CFile, "y;C;");

                for(RunVFLine){
                        for(RunVFCol){
                                for(RunCAinVFLine){
                                        for(RunCAinVFCol){

                                                if(CaP.state != 0){
                                                                esCs += CaP.fS*CaP.Cs;
                                                                eeCe += CaP.fE*Ceut;
                                                                edCd += CaP.Cd*(1.0 - CaP.fS - CaP.fE);
                                                }
                                                else{
                                                                elCl += VfP.CL;
                                                }

                                        }
                                }
                                CVF[ia] = (esCs + eeCe + edCd + elCl)/mmiy/mmix;
                                esCs = 0.0;
                                eeCe = 0.0;
                                edCd = 0.0;
                                elCl = 0.0;
                                fprintf(CFile, "\n%f;%f", Ly - dyVF*(ia + 0.5), CVF[ia] );
                        }


                }

                fclose(CFile);
                */

}

double AverageCsCdCl(double CVF[mmacrox*mmacroy][5]){

                double CAverage;
                double sumCd = 0.0, sumflCl = 0.0, sumCs = 0.0, sumCe = 0.0;
                double sumfd = 0.0, sumfs = 0.0, sumfl = 0.0, sumfE = 0.0;
                double Ce = (Teut - Tf)/ml;

                for(RunVFLine){
                for(RunVFCol){
                        CVF[ia*mmacrox + ja][0] = 0.0;
                        CVF[ia*mmacrox + ja][1] = 0.0;
                        CVF[ia*mmacrox + ja][2] = 0.0;
                        CVF[ia*mmacrox + ja][3] = 0.0;
                        CVF[ia*mmacrox + ja][4] = 0.0;
                for(RunCAinVFLine){
                for(RunCAinVFCol){

                                                if(fabs(CA[ii*mx + ji].state) >= 2){
                                                                sumfs += CaP.fS/mmix/mmiy;
                                                                sumfE += CaP.fE/mmix/mmiy;
                                                                sumfd += (1.0 - CaP.fS - CaP.fE)/mmix/mmiy;
                                                                if( CA[ii*mx + ji].fS != 0.0)
                                                                                sumCs += CA[ii*mx + ji].Cs*CA[ii*mx + ji].fS/mmiy/mmix;
                                                                if( CaP.fE != 0.0)
                                                                                sumCe += CaP.fE*Ce/mmiy/mmix;

                                                                sumCd += (CA[ii*mx + ji].Cd)*(1.0 - CaP.fS - CaP.fE)/mmix/mmiy;
                                                }
                                                else{
                                                                sumflCl += VF[ia*mmax + ja].CL/mmiy/mmix;
                                                                sumfl += 1.0/mmix/mmiy;
                                                }

                }
                }

                if(sumfl == 0.0){
                                sumfl = fLVFmin;
                                sumflCl = fLVFmin*VfP.CL;

                }


                CAverage += (sumCd + sumflCl + sumCs + sumCe);
                CVF[ia*mmacrox + ja][0] = (sumCd + sumflCl+ sumCs + sumCe);
                CVF[ia*mmacrox + ja][1] = sumCd;
                CVF[ia*mmacrox + ja][2] = sumflCl;
                CVF[ia*mmacrox + ja][3] = sumCs;
                CVF[ia*mmacrox + ja][4] = sumCe;
                sumfs = 0.0; sumfE = 0.0; sumfd = 0.0; sumCs = 0.0; sumCe = 0.0; sumCd = 0.0; sumflCl = 0.0; sumfl = 0.0;
                }
                }

                //sumCd = sumCd/(1.0 - sumfl);
                //sumCs = sumCs/(1.0 - sumfl);



                CAverage = CAverage/Lx/Ly*dxVF*dyVF;
                if(fabs(CAverage - Co) > 1.0e-10)
                        ii=ii;

                return CAverage;


}

void InitiateColors(float red[nclassestheta], float green[nclassestheta], float blue[nclassestheta]){ // Color for each class of CA){

        red[0] = 0.8;
        green[0] = 0.1;
        blue[0] = 0.6;

        red[1] = 0.2;
        green[1] = 0.1;
        blue[1] = 0.6;

        red[2] = 0.2;
        green[2] = 0.1;
        blue[2] = 0.8;

        red[3] = 0.2;
        green[3] = 0.3;
        blue[3] = 0.2;

        red[4] = 0.2;
        green[4] = 0.3;
        blue[4] = 0.4;

        red[5] = 0.2;
        green[5] = 0.3;
        blue[5] = 0.6;

        red[6] = 0.2;
        green[6] = 0.3;
        blue[6] = 0.8;

        red[7] = 0.2;
        green[7] = 0.5;
        blue[7] = 0.2;

        red[8] = 0.2;
        green[8] = 0.5;
        blue[8] = 0.4;

        red[9] = 0.2;
        green[9] = 0.5;
        blue[9] = 0.6;

        red[10] = 0.2;
        green[10] = 0.5;
        blue[10] = 0.8;

        red[11] = 0.2;
        green[11] = 0.7;
        blue[11] = 0.2;

        red[12] = 0.2;
        green[12] = 0.7;
        blue[12] = 0.4;

        red[13] = 0.2;
        green[13] = 0.7;
        blue[13] = 0.6;

        red[14] = 0.2;
        green[14] = 0.7;
        blue[14] = 0.8;

        red[15] = 0.4;
        green[15] = 0.1;
        blue[15] = 0.2;

        red[16] = 0.4;
        green[16] = 0.1;
        blue[16] = 0.4;

        red[17] = 0.4;
        green[17] = 0.1;
        blue[17] = 0.6;

        red[18] = 0.4;
        green[18] = 0.1;
        blue[18] = 0.8;

        red[19] = 0.4;
        green[19] = 0.3;
        blue[19] = 0.2;

        red[20] = 0.4;
        green[20] = 0.3;
        blue[20] = 0.4;

        red[21] = 0.4;
        green[21] = 0.3;
        blue[21] = 0.6;

        red[22] = 0.4;
        green[22] = 0.3;
        blue[22] = 0.8;

        red[23] = 0.4;
        green[23] = 0.5;
        blue[23] = 0.2;

        red[24] = 0.4;
        green[24] = 0.5;
        blue[24] = 0.4;

        red[25] = 0.4;
        green[25] = 0.5;
        blue[25] = 0.6;

        red[26] = 0.4;
        green[26] = 0.5;
        blue[26] = 0.8;

        red[27] = 0.4;
        green[27] = 0.7;
        blue[27] = 0.2;

        red[28] = 0.4;
        green[28] = 0.7;
        blue[28] = 0.4;

        red[29] = 0.4;
        green[29] = 0.7;
        blue[29] = 0.6;

        red[30] = 0.4;
        green[30] = 0.7;
        blue[30] = 0.8;

        red[31] = 0.6;
        green[31] = 0.1;
        blue[31] = 0.2;

        red[32] = 0.6;
        green[32] = 0.1;
        blue[32] = 0.4;

        red[33] = 0.6;
        green[33] = 0.1;
        blue[33] = 0.6;

        red[34] = 0.6;
        green[34] = 0.1;
        blue[34] = 0.8;

        red[35] = 0.6;
        green[35] = 0.3;
        blue[35] = 0.2;

        red[36] = 0.6;
        green[36] = 0.3;
        blue[36] = 0.4;

        red[37] = 0.6;
        green[37] = 0.3;
        blue[37] = 0.6;

        red[38] = 0.6;
        green[38] = 0.3;
        blue[38] = 0.8;

        red[39] = 0.6;
        green[39] = 0.5;
        blue[39] = 0.2;

        red[40] = 0.6;
        green[40] = 0.5;
        blue[40] = 0.4;

        red[41] = 0.6;
        green[41] = 0.5;
        blue[41] = 0.6;

        red[42] = 0.6;
        green[42] = 0.5;
        blue[42] = 0.8;

        red[43] = 0.6;
        green[43] = 0.7;
        blue[43] = 0.2;

        red[44] = 0.6;
        green[44] = 0.7;
        blue[44] = 0.4;

        red[45] = 0.6;
        green[45] = 0.7;
        blue[45] = 0.6;

        red[46] = 0.6;
        green[46] = 0.7;
        blue[46] = 0.8;

        red[47] = 0.8;
        green[47] = 0.1;
        blue[47] = 0.2;

        red[48] = 0.8;
        green[48] = 0.1;
        blue[48] = 0.4;
}





//FILE *Test = fopen( _tchdir/ "resultshi.txt", "w");






#if PRINTGRAIN
FILE *gFile;
#endif


struct CADEFiles{
                FILE *Aa;
                FILE *La;
                FILE *T;
                FILE *eg;
                FILE *es;
                FILE *CL;
                FILE *Se;
                FILE *Cd;
                FILE *dfS;

};



struct CADEFiles CADEFile;





#ifdef THERMALPROFILE
void PrintTMacro(){

        ja = 0;
        fprintf(TMacro, "\n%f;", t);


        for(ia=0; ia < mmay; ia+=1){
                fprintf(TMacro, "%f;", VF[ia*mmax + ja].T - 273.0);
        }

}

#endif




#ifdef SOLIDFRACTIONPROFILE
void PrintfSMacro(){


        fprintf(fSMacro, "\n%f;", t);

        ja = 0;
        for(ia=0; ia < mmay; ia+=1){
                fprintf(fSMacro, "%f;", VF[ia*mmax + ja].fS);
        }

}

#endif

#ifdef Eutetic
#ifdef EUTETICFRACTIONPROFILE
void PrintfEMacro(){


        fprintf(fEMacro, "\n%f;", t);

        ja = 0;
        for(ia=0; ia < mmay; ia+=1){
                fprintf(fEMacro, "%f;", VF[ia*mmax + ja].fE);
        }

}
#endif
#ifdef SUMSOLIDEUTETICFRACTIONPROFILE
void PrintfSfEMacro(){


        fprintf(fSfEMacro, "\n%f;", t);

        ja = 0;
        for(ia=0; ia < mmay; ia+=1){
                fprintf(fSfEMacro, "%f;", VF[ia*mmax + ja].fE + VF[ia*mmax + ja].fS);
        }

}
#endif
#endif

#ifdef PRINTLIQUIDFRACTIONPROFILE
void PrintfLMacro(){


        fprintf(fLMacro, "\n%f;", t);

        ja = 0;
        for(ia=0; ia < mmay; ia+=1){
                fprintf(fLMacro, "%f;", VF[ia*mmax + ja].fL);
        }

        fprintf(fLCAMacro, "\n%f;", t);

        ja = 0;
        for(ia=0; ia < mmay; ia+=1){
                fprintf(fLCAMacro, "%f;", VF[ia*mmax + ja].fLCA);
        }

}
#endif

#ifdef PRINTSOLUTEPROFILE

/*
void PrintSoluteMacro(){

double CdAverage = 0;
double CdAverageinterface;
double sumfd = 0;
double sumfdinterface = 0;
double CAverage = Co;

for(ii = 0; ii < my ; ii++){
                for(ji = 0 ; ji < mx ; ji++){
                                if(fabs(CA[ii*mx + ji].state) >= 2){
                                                CdAverage += (1.0 - CA[ii*mx + ji].fS)*CA[ii*mx + ji].Cd;
                                                sumfd += (1.0 - CA[ii*mx + ji].fS);
                                }
                                if(CA[ii*mx + ji].state >= 2){
                                                CdAverageinterface += (1.0 - CA[ii*mx + ji].fS)*CA[ii*mx + ji].Cd;
                                                sumfdinterface += (1.0 - CA[ii*mx + ji].fS);
                                }
                }
}




                if(sumfd != 0.0)
                                CdAverage = CdAverage/sumfd;

                if(sumfdinterface != 0.0)
                                CdAverageinterface = CdAverageinterface/sumfdinterface;

                CAverage = AverageCsCdCl();

                fprintf(flClFile, "\n%f;%f;%f;%f;%f;%f;%f;%f;", t, VF[0].fL, VF[ia*mmax + ja].CL, CdAverage, CdAverageinterface, VF[ia*mmax + ja].T, VF[ia*mmax + ja].T - Tf - ml*CdAverage, CAverage);




}
*/

void PrintSoluteMacro(){

                fprintf(flClFile, "\n%f;", t);

        ja = 0;
        for(ia=0; ia < mmay; ia+=1){
                fprintf(flClFile, "%f;", VF[ia*mmax + ja].CL);
        }

}

#endif

#ifdef THERMALHOPPROFILE
void PrintTMacroHop(){


        fprintf(TMacroHop, "\n%f;", t);

        ja = 0;
        for(ia= (int) ((int) (mmay/20)/2); ia < mmay; ia+= (int) (mmay/20)){
                fprintf(TMacroHop, "%f;", VF[ia*mmax + ja].T);
        }

}
#endif

#ifdef SOLIDFRACTIONHOPPROFILE
void PrintfSMacroHop(){


        fprintf(fSMacroHop, "\n%f;", t);

        ja = 0;
        for(ia= (int) ((int) (mmay/20)/2); ia < mmay; ia+= (int) (mmay/20)){
                fprintf(fSMacroHop, "%f;", VF[ia*mmax + ja].fS);
        }

}
#endif

#ifdef PRINTAverageVelocities
void PrintAverageVelocities(){
           

        fprintf(VLMacro, "\n%f;", t);
        fprintf(VdMacro, "\n%f;", t);
        fprintf(VsMacro, "\n%f;", t);
        fprintf(elVLMacro, "\n%f;", t);

        ja = 0;
        for(ia=0; ia < mmay; ia+=1){
                fprintf(VLMacro, "%.12e;", VF[ia*mmax + ja].VL);
                fprintf(VdMacro, "%.12e;", VfP.Vd);
                fprintf(VsMacro, "%.12e;", VF[ia*mmax + ja].Vs);
                //fprintf(elVLMacro, "%.12e;%.12e;", (double) elvlyn[ia], VF[ia*mmax + ja].fL*VF[ia*mmax + ja].VL);
        }
        fprintf(elVLMacro, "0.0;");
}
#endif



#ifdef PRINTGRAINS

void PrintGrains(){

                // 3.5) Variáveis para visualização no MATLAB
        float red[nclassestheta], green[nclassestheta], blue[nclassestheta]; // Color for each class of CA
        InitiateColors(red,green,blue);

        double P1[2], P2[2], P3[2], P4[2];      /* P1 = (jo , m - io) , P2 = (jo, m - jo - dm) , P3 = (j + dm , m - i - dm) , P4 = (j + dm , m - i) */
                                                /* Note que jo e io são os indices que iniciam o retângulo sólido, dm é o tamanho de cada célula, e m o tamanho da malha./*
                                                /* Houve a inversão de i e j, pois i = linha = y no MATLAB, e j = coluna = x no MATLAB */
                                                /* O desenho é feito no sentido anti-horário, iniciando no ponto P1 que é o ponto superior esquerdo*/
                                                /* O "m - i" é utilizado para a visualização ser equivalente a visualização no DOS */
        bool draw=0;



                // 3.5.1) Padronização de cores em função da orientação



                printf("\nGrao Inteiro Impresso, t = %f", t);

                draw = 0;

                fprintf (gFile, "\n\n\nt = %f", t);
                fprintf (gFile, "\n\nset (gca,'PlotBoxAspectRatio',[1 1 1],'Xlim',[0,%f],'Ylim',[0,%f])", Lx, Ly);
                fprintf (gFile, "\nhold on");
                fprintf (gFile, "\ntitle('Malha(%dX%d) dt(0.0001 s) tempo(%f s)')", mx, my, t);
                fprintf (gFile, "\nhold on");
                fprintf (gFile, "\nylabel('metros')");
                fprintf (gFile, "\nhold on");
                fprintf (gFile, "\nxlabel('metros')");
                fprintf (gFile, "\nhold on");


                for(ii=0;ii<=my-1;ii++){
                        draw = 0;
                        for(ji=0;ji<=mx-1;ji++){
                                if( draw == 0 && ji != 0 && CA[ii*mx + ji].state!=0 && fabs(CA[ii*mx + ji].state) != fabs(CA[ii*mx + ji - 1].state) ){
                                                P1[0] = ji;
                                                P1[1] = my - ii;
                                                P2[0] = ji;
                                                P2[1] = my - ii - 1;

                                                draw=1;
                                }

                                if( draw == 0 && ji == 0 && CA[ii*mx + ji].state!=0){
                                                P1[0] = ji;
                                                P1[1] = my - ii;
                                                P2[0] = ji;
                                                P2[1] = my - ii - 1;

                                                draw=1;
                                }


                                if( draw == 1 && ji != mx - 1 && fabs(CA[ii*mx + ji].state) != fabs(CA[ii*mx + ji+ 1].state) ){
                                                P3[0] = ji + 1;
                                                P3[1] = my - ii - 1;
                                                P4[0] = ji + 1;
                                                P4[1] = my - ii;

                                                draw=0;

                                                #ifndef PRINTGrainsClusters
                                                fprintf(gFile, "\nfill([%f %f %f %f],[%f %f %f %f],[%.1f %.1f %.1f],'LineStyle','none')", P1[0]*dxCA, P2[0]*dxCA, P3[0]*dxCA, P4[0]*dxCA, P1[1]*dxCA, P2[1]*dxCA, P3[1]*dxCA, P4[1]*dxCA, red[ clusters[(int) fabs(CA[ii*mx + ji].state) - 2].thetaclass], green[clusters[(int) fabs(CA[ii*mx + ji].state) - 2].thetaclass], blue[clusters[(int) fabs(CA[ii*mx + ji].state) - 2].thetaclass]);
                                                #endif
                                                #ifdef PRINTGrainsClusters
                                                fprintf(gFile, "\nfill([%f %f %f %f],[%f %f %f %f],[%.1f %.1f %.1f],'LineStyle','none')", P1[0]*dxCA, P2[0]*dxCA, P3[0]*dxCA, P4[0]*dxCA, P1[1]*dxCA, P2[1]*dxCA, P3[1]*dxCA, P4[1]*dxCA, red[ clusters[ClusterP].thetaclass], green[clusters[ClusterP].thetaclass], blue[clusters[ClusterP].thetaclass]);
                                                #endif
                                }

                                if( draw == 1 && ji == mx - 1){
                                                P3[0] = ji + 1;
                                                P3[1] = my - ii - 1;
                                                P4[0] = ji + 1;
                                                P4[1] = my - ii;

                                                draw=0;

                                                #ifndef PRINTGrainsClusters
                                                fprintf(gFile, "\nfill([%f %f %f %f],[%f %f %f %f],[%.1f %.1f %.1f],'LineStyle','none')", P1[0]*dxCA, P2[0]*dxCA, P3[0]*dxCA, P4[0]*dxCA, P1[1]*dxCA, P2[1]*dxCA, P3[1]*dxCA, P4[1]*dxCA, red[ clusters[(int) fabs(CA[ii*mx + ji].state) - 2].thetaclass], green[clusters[(int) fabs(CA[ii*mx + ji].state) - 2].thetaclass], blue[clusters[(int) fabs(CA[ii*mx + ji].state) - 2].thetaclass]);
                                                #endif
                                                #ifdef PRINTGrainsClusters
                                                fprintf(gFile, "\nfill([%f %f %f %f],[%f %f %f %f],[%.1f %.1f %.1f],'LineStyle','none')", P1[0]*dxCA, P2[0]*dxCA, P3[0]*dxCA, P4[0]*dxCA, P1[1]*dxCA, P2[1]*dxCA, P3[1]*dxCA, P4[1]*dxCA, red[ clusters[ClusterP].thetaclass], green[clusters[ClusterP].thetaclass], blue[clusters[ClusterP].thetaclass]);
                                                #endif
                                }
                        }
                }


}

#endif

#ifdef PRINTCONTORNO


void PrintContorno(){

                float red[nclassestheta], green[nclassestheta], blue[nclassestheta]; // Color for each class of CA
        InitiateColors(red,green,blue);
        int P1[2], P2[2], P3[2], P4[2];

       printf("\nContorno Impresso, t = %f\n", t);


                        fprintf (cFile, "\n\n\nt = %f", t);
                        fprintf (cFile, "\n\nset (gca,'PlotBoxAspectRatio',[1 1 1],'Xlim',[0,%f],'Ylim',[0,%f])", Lx, Ly);
                        fprintf (cFile, "\nhold on");
                        fprintf (cFile, "\ntitle('Malha(%dX%d) dt(0.001 s) tempo(%f s)')", mx, my, t);
                        fprintf (cFile, "\nhold on");
                        fprintf (cFile, "\nylabel('metros')");
                        fprintf (cFile, "\nhold on");
                        fprintf (cFile, "\nxlabel('metros')");
                        fprintf (cFile, "\nhold on");


                        //fprintf (centroFile, "\n\n\nt = %f", t);


                        for(i=1;i<my-1;i++){
                                for(j=1;j<=mx-1;j++){
                                        if( CA[i*mx + j].state>=1){
                                                P1[0] = j;
                                                P1[1] = my - i;
                                                P2[0] = j;
                                                P2[1] = my - i - 1;

                                                P3[0] = j + 1;
                                                P3[1] = my - i - 1;
                                                P4[0] = j + 1;
                                                P4[1] = my - i;

                                                // printf("\n L[%d,%d] = %e", i, j, volmicro[i*mx + j].L);


                                                fprintf(cFile, "\nfill([%f %f %f %f],[%f %f %f %f],[%f %f %f],'LineStyle','none')", P1[0]*dxCA, P2[0]*dxCA, P3[0]*dxCA, P4[0]*dxCA, P1[1]*dyCA, P2[1]*dyCA, P3[1]*dyCA, P4[1]*dyCA, red[clusters[(int) fabs(CA[i*mx + j].state) - 1].thetaclass]*0.1, green[clusters[(int) fabs(CA[i*mx + j].state) - 1].thetaclass]*0.1, blue[clusters[(int) fabs(CA[i*mx + j].state) - 1].thetaclass]*0.1);


                                        }

                                        /*
                                        if( volmicro[i*mx + j].solidfraction == 1.0){

                                                // Cruz para indicação do centro de crescimento

                                                P1[0] = j - 1;
                                                P1[1] = my - i + 1;
                                                P2[0] = j;
                                                P2[1] = my - i + 1;

                                                P3[0] = j + 2;
                                                P3[1] = my - i - 2;
                                                P4[0] = j + 1;
                                                P4[1] = my - i - 2;

                                                fprintf(centroFile, "\nfill([%f %f %f %f],[%f %f %f %f],[%f %f %f],'LineStyle','none')", P1[0]*dxCA, P2[0]*dxCA, P3[0]*dxCA, P4[0]*dxCA, P1[1]*dyCA, P2[1]*dyCA, P3[1]*dyCA, P4[1]*dyCA, red[registro[(int) fabs(volmicro[i*mx + j].state) - 1].classe]*0.1, green[registro[(int) fabs(volmicro[i*mx + j].state) - 1].classe]*0.1, blue[registro[(int) fabs(volmicro[i*mx + j].state) - 1].classe]*0.1);

                                                P1[0] = j + 1;
                                                P1[1] = my - i + 1;
                                                P2[0] = j + 2;
                                                P2[1] = my - i + 1;

                                                P3[0] = j;
                                                P3[1] = my - i - 2;
                                                P4[0] = j - 1;
                                                P4[1] = my - i - 2;

                                                fprintf(centroFile, "\nfill([%f %f %f %f],[%f %f %f %f],[%f %f %f],'LineStyle','none')", P1[0]*dxCA, P2[0]*dxCA, P3[0]*dxCA, P4[0]*dxCA, P1[1]*dyCA, P2[1]*dyCA, P3[1]*dyCA, P4[1]*dyCA, red[registro[(int) fabs(volmicro[i*mx + j].state) - 1].classe]*0.1, green[registro[(int) fabs(volmicro[i*mx + j].state) - 1].classe]*0.1, blue[registro[(int) fabs(volmicro[i*mx + j].state) - 1].classe]*0.1);
                                        }
                                        */
                                }

                        }
}

#endif



double CalculateSe( double La ){

        return (4.0/M_PI)*La;

}



double Calculatefg(double Aa){


        return Aa;


}


void CountLaAa(double *La, double *Aa){

        *La = 0.0;
        *Aa = 0.0;


        struct MicroMatrix *CellPointer ;

        for(RunCAinVFLine){
                for (RunCAinVFCol){
                                if( t > 11.576 )
                                                ii=ii;
                        CellPointer = &CA[ii*mx + ji];
                        if( fabs(CA[ii*mx + ji].state)  >= 2.0 ){

                                *Aa = *Aa + 1;
                                ii=ii;


                                if( ii != 0.0 ){
                                        if( fabs((CellPointer - mx)->state) == 0.0){
                                                *La = *La + 1;
                                        }

                                }

                                if( ii != my - 1 ){
                                        if( fabs((CellPointer + mx)->state) == 0.0){
                                                *La = *La + 1;
                                        }

                                }

                                if( ji != 0.0 ){
                                        if( fabs((CellPointer - 1)->state) == 0.0){
                                                *La = *La + 1;
                                        }

                                }

                                if( ji != mx - 1 ){
                                        if( fabs((CellPointer + 1)->state) == 0.0){
                                               *La = *La + 1;
                                        }

                                }

                                if( ii == 0 || ii == my - 1 || ji == 0 || ji == mx - 1)
                                                ii=ii;


                        }




                }
        }

        if( *La != 0.0 ){
                ii =ii;
        }
        *La = *La/mmicrox/mmicroy/dxCA;
        *Aa = *Aa*dxCA*dyCA/(dxVF*dyVF);

}




void InitiatePrintCADE(){

                CADEFile.Aa = fopen("Aa.dat", "w");
                CADEFile.La = fopen("La.dat", "w");
                CADEFile.T = fopen("T.dat", "w");
                CADEFile.eg = fopen("eg.dat", "w");
                CADEFile.es = fopen("es.dat", "w");
                CADEFile.CL = fopen("CL.dat", "w");
                CADEFile.Se = fopen("Se.dat", "w");
                CADEFile.Cd = fopen("Cd.dat", "w");
                CADEFile.dfS = fopen("dfS.dat","w");

                //fprintf(CADEFile.Aa, "t;Aa");
                fprintf(CADEFile.La, "t;La");
                fprintf(CADEFile.T, "t;T");
                fprintf(CADEFile.eg, "t;eg");
                fprintf(CADEFile.es, "t;fS");
                fprintf(CADEFile.CL, "t;CL");
                fprintf(CADEFile.Se, "t;Se");
                fprintf(CADEFile.Cd, "t;CdAv;CdIntAve;Cdmin;Cdmax;Cdeq;Co");
                fprintf(CADEFile.dfS, "t;dfS");





}

void AverageCdforVF(){

                Cdmin = 100.0;
                Cdmax = 0.0;
                double AvCd = 0.0;
                double AvCdInterface = 0.0;
                double CdCACounter = 0;
                double CdCAInterfaceCounter = 0;
                VFP = &VfP;

                double AvCs = 0.0;
                double CsCACounter = 0.0;

                for(RunCAinVFLine){
                                for(RunCAinVFCol){
                                                CAP = &CA[ii*mx + ji];
                                                if( fabs(CAP->state) >= 2.0){
                                                                CdCACounter += (1.0 - CAP->fS);
                                                                CsCACounter += CAP->fS;
                                                                AvCd += (1.0 - CAP->fS)*CAP->Cd;
                                                                AvCs += CAP->fS*CAP->Cs;
                                                                if( CAP->Cd < Cdmin){
                                                                                Cdmin = CAP->Cd;
                                                                }
                                                                if( CAP->Cd > Cdmax){
                                                                                Cdmax = CAP->Cd;
                                                                }


                                                }

                                                if(CAP->state >= 2 ){
                                                                CdCAInterfaceCounter += (1.0 - CAP->fS);
                                                                AvCdInterface += (1.0 - CAP->fS)*CAP->Cd;

                                                }
                                }
                }

                if(CdCACounter !=0){
                                CdAve = AvCd/CdCACounter;
                }
                else{
                                CdAve = Co;
                                Cdmin = Co;
                                Cdmax = Co;
                                CdInterfaceAve = Co;
                }
                
                if(CsCACounter != 0.0){
                                CsAve = AvCs/CsCACounter;
                }
                else{
                                CsAve = 0.0;
                }
                if(CdCAInterfaceCounter!=0)
                                CdInterfaceAve = AvCdInterface/CdCAInterfaceCounter;
                else
                                CdInterfaceAve = Co;

                if( CdInterfaceAve > Cdmax)
                                ia=ia;

}

void PrintCADE(){

                double La, Aa  ;
                ia = 0;
                ja = 0;
                VFP = &VfP;

                CountLaAa(&La,&Aa);
                AverageCdforVF();

                fprintf(CADEFile.Aa, "\n%f;%e;", t, Aa);
                fprintf(CADEFile.La, "\n%f;%e;", t, La);
                fprintf(CADEFile.T, "\n%f;%e;", t, VFP->T);
                fprintf(CADEFile.eg, "\n%f;%e;", t, 1.0 - VFP->fL);
                fprintf(CADEFile.es, "\n%f;%e;", t, VFP->fS);
                fprintf(CADEFile.CL, "\n%f;%e;", t, VFP->CL);
                
                fprintf(CADEFile.Se, "\n%f;%e;", t, /*CalculateSe(La)*/ JDdlVF/(CdAve - VfP.CL)*(del_dl/DiffLiq)/dt );
                fprintf(CADEFile.Cd, "\n%f;%e;%e;%e;%e;%e" , t, Tf + ml*CdAve, Tf + ml*CdInterfaceAve, Tf + ml*Cdmin , Tf + ml*Cdmax, VF[ia*mmax + ja].T);
                fprintf(CADEFile.Cd, ";%e", (VF[ia*mmax + ja].fS*CsAve + (1.0 - VF[ia*mmax + ja].fL - VF[ia*mmax + ja].fS)*CdAve + VF[ia*mmax + ja].fL*VF[ia*mmax + ja].CL)*ml + Tf);
                

}





void PrintCADEComparisonResults(){

                //fprintf();
}

void PrintProfiles(){

#ifdef THERMALPROFILE
PrintTMacro();
#endif

#ifdef SOLIDFRACTIONPROFILE
PrintfSMacro();
#endif

#ifdef Eutetic
#ifdef EUTETICFRACTIONPROFILE
PrintfEMacro();
#endif
#ifdef SUMSOLIDEUTETICFRACTIONPROFILE
PrintfSfEMacro();
#endif
#endif

#ifdef PRINTSOLUTEPROFILE
PrintSoluteMacro();
#endif

#ifdef PRINTLIQUIDFRACTIONPROFILE
PrintfLMacro();
#endif

#ifdef THERMALHOPPROFILE
PrintTMacroHop();
#endif

#ifdef SOLIDFRACTIONHOPPROFILE
PrintfSMacroHop();
#endif

#ifdef PRINTAverageVelocities
PrintAverageVelocities();
#endif

#ifdef PRINTSeProfile
fprintf(Evolution_Se, "\n%f;", t);
fprintf(Evolution_Ne, "\n%f;", t);
fprintf(Evolution_Nuc, "\n%f;", t);
for(RunVFLine){
        for(RunVFCol){
                fprintf(Evolution_Se, "%f;", VfP.Se);
                fprintf(Evolution_Ne, "%f;", pow(VfP.Ne/dyVF/dxVF, 3.0/2.0)*sqrt(M_PI/6.0));
                fprintf(Evolution_Nuc, "%d;", VfP.Nuc);
        }
}
#endif

}


void InitialVisual(){

        FILE *StartFile;
        FILE *dTFile;

        StartFile = fopen ("inicializacao.txt", "w");
        dTFile = fopen ("substratos.txt", "w");

        int P1[2], P2[2], P3[2], P4[2];
        int contadordraw = 0;
                printf("\nCarregando Visualizacao de inicializacao...");

                 // Visualização no MATLAB

                        fprintf (StartFile, "\n\n\ncontadorS = %d", contadorS);
                        fprintf (StartFile, "\n\n\ncontadorV = %d", contadorV);
                        fprintf (StartFile, "\n\n\nt = %d", 8);
                        fprintf (StartFile, "\n\nset (gca,'PlotBoxAspectRatio',[1 1 1],'Xlim',[0,%f],'Ylim',[0,%f])", Lx, Ly);
                        fprintf (StartFile, "\nhold on");
                        fprintf (StartFile, "\ntitle('Mesh(%dX%d) dt(%f s) t(%f s)')", mx, my, dt,0.01);
                        fprintf (StartFile, "\nhold on");
                        fprintf (StartFile, "\nylabel('metros')");
                        fprintf (StartFile, "\nhold on");
                        fprintf (StartFile, "\nxlabel('metros')");
                        fprintf (StartFile, "\nhold on");


                        for(ii=0;ii<my;ii++){
                                for(ji=0;ji<mx;ji++){
                                        if( CA[ii*mx + ji].dT>=0.0){

                                                fprintf(dTFile, "\n%f;", CA[ii*mx+ji].dT);

                                                P1[0] = ji;
                                                P1[1] = my - ii;
                                                P2[0] = ji;
                                                P2[1] = my - ii - 1;

                                                P3[0] = ji + 1;
                                                P3[1] = my - ii - 1;
                                                P4[0] = ji + 1;
                                                P4[1] = my - ii;


                                                fprintf(StartFile, "\nfill([%f %f %f %f],[%f %f %f %f],[0.1 0.4 0.8],'LineStyle','none')", P1[0]*dxCA, P2[0]*dxCA, P3[0]*dxCA, P4[0]*dxCA, P1[1]*dxCA, P2[1]*dxCA, P3[1]*dxCA, P4[1]*dxCA);
                                                contadordraw++  ;
                                        }
                                }
                        }
                        fprintf (StartFile, "\n\n\ncontadordraw = %d", contadordraw);


                        fclose(StartFile);

}


void InitiateFiles(){
#ifdef THERMALPROFILE
        TMacro = fopen("Evolution_T.dat", "w");
        fprintf(TMacro, "t;");
        for(ia=0; ia<mmay; ia++){
                fprintf(TMacro, "y = %.3f m;", Ly - (ia + 0.5)*dyVF);
        }
        #endif

        #ifdef SOLIDFRACTIONPROFILE
        fSMacro = fopen("Evolution_fS.dat","w");
        fprintf(fSMacro, "t;");
        for(ia=0; ia<mmay; ia++){
                fprintf(fSMacro, "y = %.3f m;", Ly - (ia + 0.5)*dyVF);
        }
        #endif

        #ifdef Eutetic
        #ifdef EUTETICFRACTIONPROFILE
        fEMacro = fopen("Evolution_fE.dat","w");
        fprintf(fEMacro, "t;");
        for(ia=0; ia<mmay; ia++){
                fprintf(fEMacro, "y = %.3f m;", Ly - (ia + 0.5)*dyVF);
        }
        #endif
        #ifdef SUMSOLIDEUTETICFRACTIONPROFILE
        fSfEMacro = fopen("Evolution_fSfE.dat","w");
        fprintf(fSfEMacro, "t;");
        for(ia=0; ia<mmay; ia++){
                fprintf(fSfEMacro, "y = %.3f m;", Ly - (ia + 0.5)*dyVF);
        }
        #endif
        #endif

        #ifdef THERMALHOPPROFILE
        TMacroHop = fopen("Evolution_T_Hop.dat","w");
        fprintf(TMacroHop, "t;");
        for(ia= (int) ((int) (mmay/20)/2); ia < mmay; ia+= max( 1, (int) (mmay/20) ) ) {
                fprintf(TMacroHop, "y = %.3f m;", Ly - (ia + 0.5)*dyVF);
        }
        #endif

        #ifdef PRINTSOLUTEPROFILE
        flClFile = fopen("Evolution_CL.dat","w");
        fprintf(flClFile, "t;");
        for(ia=0; ia<mmay; ia++){
                fprintf(flClFile, "y = %.3f m;", Ly - (ia + 0.5)*dyVF);
        }

        #endif

        #ifdef PRINTLIQUIDFRACTIONPROFILE
        fLMacro = fopen("Evolution_fL.dat","w");
        fprintf(fLMacro, "t;");
        for(ia=0; ia<mmay; ia++){
                fprintf(fLMacro, "y = %.3f m;", Ly - (ia + 0.5)*dyVF);
        }


        fLCAMacro = fopen("Evolution_fLCA.dat","w");
        fprintf(fLCAMacro, "t;");
        for(ia=0; ia<mmay; ia++){
                fprintf(fLCAMacro, "y = %.3f m;", Ly - (ia + 0.5)*dyVF);
        }
        #endif

        #ifdef PRINTAverageVelocities
        VLMacro = fopen("Evolution_VL.dat", "w");
        fprintf(VLMacro, "t;");
        for(ia=0; ia<mmay; ia++){
                fprintf(VLMacro, "y = %.3f m;", Ly - (ia + 0.5)*dyVF);
        }

        VdMacro = fopen("Evolution_Vd.dat", "w");
        fprintf(VdMacro, "t;");
        for(ia=0; ia<mmay; ia++){
                fprintf(VdMacro, "y = %.3f m;", Ly - (ia + 0.5)*dyVF);
        }

        VsMacro = fopen("Evolution_Vs.dat", "w");
        fprintf(VsMacro, "t;");
        for(ia=0; ia<mmay; ia++){
                fprintf(VsMacro, "y = %.3f m;", Ly - (ia + 0.5)*dyVF);
        }

        elVLMacro = fopen("Evolution_eLVL.dat", "w");
        fprintf(elVLMacro, "t;");
        for(ia=0; ia<mmay; ia++){
                fprintf(elVLMacro, "%f;%f;", Ly - ia*dyVF, Ly - dyVF/2.0 - ia*dyVF);
        }
        fprintf(elVLMacro, "0.0;");

        #endif

        #ifdef PRINTGRAINS
        gFile = fopen("Grains.txt", "w");
        #endif

        #ifdef PRINTCONTORNO
        cFile = fopen("Contorno.txt", "w");
        #endif

        #ifdef PRINTCADERESULTS
        InitiatePrintCADE();
        #endif

        
        
        #ifdef PRINTSURFACECD
        SurfaceCd = fopen("SurfaceCd.m", "w");
        fprintf(SurfaceCd, "X = [");
        SurfacefS = fopen("SurfacefS.m", "w");
        fprintf(SurfacefS, "X = [");
        SurfaceVd = fopen("SurfaceVd.m", "w");
        fprintf(SurfaceVd, "X = [");

        for(RunCACol){
                fprintf(SurfaceCd, "%.10e",(0.5 + ji)*dxCA);
                fprintf(SurfacefS, "%.10e",(0.5 + ji)*dxCA);
                fprintf(SurfaceVd, "%.10e",(0.5 + ji)*dxCA);
                if(ji<mx-1){
                        fprintf(SurfaceCd,",");
                        fprintf(SurfacefS,",");
                        fprintf(SurfaceVd,",");
                }
                else{
                        fprintf(SurfaceCd,"]");
                        fprintf(SurfacefS,"]");
                        fprintf(SurfaceVd,"]");
                }

        }

        fprintf(SurfaceCd, "\nY = [");
        fprintf(SurfacefS, "\nY = [");
        fprintf(SurfaceVd, "\nY = [");
        for(RunCALine){
                fprintf(SurfaceCd, "%.10e", (0.5 + ii)*dyCA);
                fprintf(SurfacefS, "%.10e", (0.5 + ii)*dyCA);
                fprintf(SurfaceVd, "%.10e", (0.5 + ii)*dyCA);
                if(ii<my-1){
                        fprintf(SurfaceCd,",");
                        fprintf(SurfacefS,",");
                        fprintf(SurfaceVd,",");

                }
                else{
                        fprintf(SurfaceCd,"]");
                        fprintf(SurfacefS,"]");
                        fprintf(SurfaceVd,"]");
                }
        }
        #endif

        efvf = fopen("efvfyn.dat", "w");
        elvl = fopen("elvlyn.dat", "w");


        fprintf(efvf, "t;");
        fprintf(elvl, "t;");

        for(RunVFLine){
                fprintf(efvf, "y = %.3f m;", Ly - (ia + 0.5)*dyVF);
                fprintf(elvl, "y = %.3f m;", Ly - (ia + 0.5)*dyVF);
        }

        #ifdef PRINTSeProfile
        Evolution_Se = fopen("Evo_Se.dat","w");
        Evolution_Ne = fopen("Evo_Ne.dat","w");
        Evolution_Nuc = fopen("Evo_Nuc.dat","w");
        fprintf(Evolution_Se, "t;");
        fprintf(Evolution_Ne, "t;");
        fprintf(Evolution_Nuc, "t;");
        for(RunVFLine){
                fprintf(Evolution_Se, "y = %.3f m;", Ly - (ia + 0.5)*dyVF);
                fprintf(Evolution_Ne, "y = %.3f m;", Ly - (ia + 0.5)*dyVF);
                fprintf(Evolution_Nuc, "y = %.3f m;", Ly - (ia + 0.5)*dyVF);
        }
        #endif



        #ifdef PRINTVERTICALVALUES
        VertCL = fopen("Vert_CL.dat", "w");
        fprintf(VertCL, "t;");
        for(RunInverseVFLine){
                fprintf(VertCL, "VF%d;", ia);
        }

        VertVL = fopen("Vert_VL.dat", "w");
        fprintf(VertVL, "t;");
        for(RunInverseVFLine){
                fprintf(VertVL, "VF%d;", ia);
        }

        VertVd = fopen("Vert_Vd.dat", "w");
        fprintf(VertVd, "t;");
        for(RunInverseVFLine){
                fprintf(VertVd, "VF%d;", ia);
        }

        VertVs = fopen("Vert_Vs.dat", "w");
        fprintf(VertVs, "t;");
        for(RunInverseVFLine){
                fprintf(VertVs, "VF%d;", ia);
        }

        VertfS = fopen("Vert_fS.dat", "w");
        fprintf(VertfS, "t;");
        for(RunInverseVFLine){
                fprintf(VertfS, "VF%d;", ia);
        }

        VertflVl = fopen("Vert_fLVL.dat", "w");
        fprintf(VertflVl, "t;");
        for(RunInverseVFLine){
                fprintf(VertflVl, "VF%d;", ia);
        }

        VertfsVs = fopen("Vert_fsVs.dat", "w");
        fprintf(VertfsVs, "t;");
        for(RunInverseVFLine){
                fprintf(VertfsVs, "VF%d;", ia);
        }

        VertfdVd = fopen("Vert_fdVd.dat", "w");
        fprintf(VertfdVd, "t;");
        for(RunInverseVFLine){
                fprintf(VertfdVd, "VF%d;", ia);
        }

        VertSe = fopen("Vert_Se.dat", "w");
        fprintf(VertSe, "t;");
        for(RunInverseVFLine){
                fprintf(VertSe, "VF%d;", ia);
        }

        VertT = fopen("Vert_T.dat", "w");
        fprintf(VertT, "t;");
        for(RunInverseVFLine){
                fprintf(VertT, "VF%d;", ia);
        }

        VertCs = fopen("Vert_Cs.dat", "w");
        fprintf(VertCs, "t;");
        for(RunInverseVFLine){
                fprintf(VertCs, "VF%d;", ia);
        }

        VertC = fopen("Vert_C.dat", "w");
        fprintf(VertC, "t;");
        for(RunInverseVFLine){
                fprintf(VertC, "VF%d;", ia);
        }
        #endif

        #ifdef PRINTCALL
        Evolution_Call = fopen("Evo_Call.dat","w");
        fprintf(Evolution_Call, "t;Call;");
        #endif



}

#ifdef PRINTVERTICALVALUES
void PrintVerticalValues(double CVF[mmacrox*mmacroy][5]){

double CsAv[mmacrox*mmacroy];
double CAv[mmacrox*mmicroy];



        ja = 0;
        fprintf(VertCL, "\n%f;", t);
        for(RunInverseVFLine){
                fprintf(VertCL, "%f;", VF[ia*mmax + ja].CL);
        }

        ja = 0;
        fprintf(VertVL, "\n%f;", t);
        for(RunInverseVFLine){
                fprintf(VertVL , "%f;", VF[ia*mmax + ja].VL);
        }

        ja = 0;
        fprintf(VertVd, "\n%f;", t);
        for(RunInverseVFLine){
                fprintf(VertVd , "%f;", VF[ia*mmax + ja].Vd);
        }

        ja = 0;
        fprintf(VertVs, "\n%f;", t);
        for(RunInverseVFLine){
                fprintf(VertVs , "%f;", VF[ia*mmax + ja].Vs);
        }

        ja = 0;
        fprintf(VertfS, "\n%f;", t);
        for(RunInverseVFLine){
                fprintf(VertfS, "%f;", VF[ia*mmax + ja].fS);
        }

        ja = 0;
        fprintf(VertflVl, "\n%f;", t);
        for(RunInverseVFLine){
                fprintf(VertflVl, "%f;", VF[ia*mmax + ja].fL*VfP.VL);
        }

        ja = 0;
        fprintf(VertfsVs, "\n%f;", t);
        for(RunInverseVFLine){
                fprintf(VertfsVs, "%f;", VF[ia*mmax + ja].fS*VfP.Vs);
        }

        ja = 0;
        fprintf(VertfdVd, "\n%f;", t);
        for(RunInverseVFLine){
                fprintf(VertfdVd, "%f;", (1.0 - VfP.fL - VF[ia*mmax + ja].fS)*VfP.Vd);
        }

        ja = 0;
        fprintf(VertSe, "\n%f;", t);
        for(RunInverseVFLine){
                fprintf(VertSe, "%f;", VF[ia*mmax + ja].Se);
        }

        ja = 0;
        fprintf(VertT, "\n%f;", t);
        for(RunInverseVFLine){
                fprintf(VertT, "%f;", VF[ia*mmax + ja].T);
        }

        ja = 0;
        fprintf(VertCs, "\n%f;", t);
        //Csolidaverage(CsAv);
        AverageCsCdCl(CVF);
        for(RunInverseVFLine){
                if(VF[ia].fS != 0.0)
                        fprintf(VertCs, "%f;", CVF[ia][3]/VF[ia].fS);
                else
                        fprintf(VertCs, "%f;", 0.0);
        }

        ja = 0;
        fprintf(VertC, "\n%f;", t);
        //Caverage(CAv);
        for(RunInverseVFLine){
                fprintf(VertC, "%f;", CVF[ia][0]);

        }




}
#endif

#ifdef PRINTCALL
void PrintCall(double CVF[mmacrox*mmacroy][5]){


        double Call;
        double Comp;


        ja = 0;
        fprintf(Evolution_Call, "\n%f;", t);
        Call = AverageCsCdCl(CVF);
        if( fabs(Call - Co)>0.0000001 )
                Comp = Call - Co;
        fprintf(Evolution_Call, "%f;", Call);
}

#endif


#ifdef PRINTSURFACECD
int indexSurfaceCd = 0;


void PrintSurfaceCd(){

                indexSurfaceCd ++;

                fprintf(SurfaceCd,"\n\n\n\n t = %f", t);
                fprintf(SurfaceCd, "\nZ%d = [", indexSurfaceCd);

                fprintf(SurfacefS,"\n\n\n\n t = %f", t);
                fprintf(SurfacefS, "\nZ%d = [", indexSurfaceCd);

                fprintf(SurfaceVd,"\n\n\n\n t = %f", t);
                fprintf(SurfaceVd, "\nZ%d = [", indexSurfaceCd);

                for(RunInverseVFLine){
                for(RunVFCol){
                for(RunInverseCAinVFLine){

                                for(RunCAinVFCol){
                                                if(CaP.state == 0.0){
                                                                //fprintf(SurfaceCd, "%.5e", VF[0].T);
                                                                //fprintf(SurfaceCd, "%.5e", VF[0].CL);
                                                                fprintf(SurfaceCd, "%.5e", VfP.CL);
                                                                fprintf(SurfacefS, "%.5e", 0.0);
                                                                fprintf(SurfaceVd, "%.5e", 0.0);
                                                }
                                                else{

                                                #ifdef HoldBoundaryCondition

                                                        if(BoundaryConditionHolder[ji].on == 0){

                                                                BoundaryConditionHolder[ji].on = 1;
                                                                if(CaP.fS >= fSmax){
                                                                        BoundaryConditionHolder[ji].edVd = 0.0;
                                                                }
                                                                else{
                                                                        //BoundaryConditionHolder[ji].edVd = (1.0 - VfP.fS - VfP.fE - VfP.fL)*(VfP.Vd - clusters[ClusterP].V[1]);
                                                                        BoundaryConditionHolder[ji].edVd = edVdVsInterpol()/(1.0 - fLInterpol());
                                                                }
                                                        }



                                #endif

                                                                //fprintf(SurfaceCd, "%.5e", Tf + ml*CA[ii*mx+ji].Cd);
                                                                //fprintf(SurfaceCd, "%.5e", CA[ii*mx+ji].Cd);
                                                                fprintf(SurfaceCd, "%.5e", CA[ii*mx+ji].Cd);
                                                                fprintf(SurfacefS, "%.5e", CA[ii*mx+ji].fS);
                                                                fprintf(SurfaceVd, "%.5e", BoundaryConditionHolder[ji].edVd/(1.0 - CaP.fS - CaP.fE));

                                                                if(CaP.state < 0.0)
                                                                                ii=ii;
                                        }

                                        if(ji<mx-1){
                                                fprintf(SurfaceCd,",");
                                                fprintf(SurfacefS,",");
                                                fprintf(SurfaceVd,",");

                                        }
                                }

                                if(ii>0){
                                                fprintf(SurfaceCd,";");
                                                fprintf(SurfacefS,";");
                                                fprintf(SurfaceVd,";");

                                }
                                else{
                                                fprintf(SurfaceCd,"]");
                                                fprintf(SurfacefS,"]");
                                                fprintf(SurfaceVd,"]");
                                }

                }
                }
                }
}

#endif









