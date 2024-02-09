#include <stdio.h>
#include <math.h>

//---------------------------------------------------------------------------
#pragma hdrstop

#include "InterpolationFunctions.h"
//---------------------------------------------------------------------------
#pragma package(smart_init)

 

// Main Includes
#include "GlobalDefinitions.h"
#include "GlobalVariables.h"
#include "SetupDatabase.h"

#include "ExtraFunctions.h"

#ifdef ThermoCouples

void SaveFrontPosition_TGrad(int save_header){
        FILE *filepointer;
        double Yposition;
        int getout;

        for(RunCALine){
                for(RunCACol){

                }
        }

        if(save_header == 1)
        {
                remove("Front.dat");
                filepointer = fopen("Front.dat", "a+");
                fprintf(filepointer, "t(s); Front Position (m); Grad T (K/m)\n");
        }
        else
                filepointer = fopen("Front.dat","a+");

        Yposition = (my - ii)*dyCA;

}

void SaveLiquidusPos(){

}

struct ThermoCoupleData ReadThermoCouple(){
        FILE *entrada;
        struct ThermoCoupleData Data;
       Data.PosInicial[0] = Pos0;
       Data.PosInicial[1] = Pos1;
       Data.PosInicial[2] = Pos2;
       Data.PosInicial[3] = Pos3;
       Data.PosInicial[4] = Pos4;
       Data.PosInicial[5] = Pos5;
       Data.PosInicial[6] = Pos6;
                     //Posi��o dos termopares

       int PosicaoCurvasResfriamento;         //Marcador da posi��o atual na curva de resfriamento

       char linha_Lixo [100];

       entrada = fopen("Cooling_Curves.dat","r");

       fgets (linha_Lixo,100,entrada); //L� toda primeira linha (cabe�alho)
       PosicaoCurvasResfriamento = 0;

       while ( fscanf (entrada,"%lf %lf %lf %lf %lf %lf %lf %lf",
                               &Data.TimeCurves[PosicaoCurvasResfriamento],
                               &Data.Curvas_Resf[0][PosicaoCurvasResfriamento],
                               &Data.Curvas_Resf[1][PosicaoCurvasResfriamento],
                               &Data.Curvas_Resf[2][PosicaoCurvasResfriamento],
                               &Data.Curvas_Resf[3][PosicaoCurvasResfriamento],
                               &Data.Curvas_Resf[4][PosicaoCurvasResfriamento],
                               &Data.Curvas_Resf[5][PosicaoCurvasResfriamento],
                               &Data.Curvas_Resf[6][PosicaoCurvasResfriamento])
              != EOF) //L� cada linha de n�mero
       PosicaoCurvasResfriamento++;
       Data.TimeCurves[PosicaoCurvasResfriamento]=-1; //Marca o final do vetor
       PosicaoCurvasResfriamento = 0;

       fclose (entrada);

       return Data;
}

void RefreshVFThermoCoupleData(double ti, struct ThermoCoupleData Data){
        double yLinha;
        double dT_dz;
        double dT_dt;
        double TTemp;
        int zLinha;
        int PositionZ;
        int PosicaoCurvasResfriamento = 0;
        while(ti > Data.TimeCurves[PosicaoCurvasResfriamento+1])
           PosicaoCurvasResfriamento++;

        if (Data.TimeCurves[PosicaoCurvasResfriamento+1] < 0)
         {
          printf ("Tentando ler mais medidas de T do que o dispon�vel");
         }

        ja = 0;
        for(RunVFLine){

                yLinha = Ly - (ia + 0.5)*dyVF;
                PositionZ = 0;
                while (yLinha > Data.PosInicial[PositionZ+1])
                PositionZ++;

                dT_dz =
                    0.5*( (Data.Curvas_Resf[PositionZ+1][PosicaoCurvasResfriamento] -
                           Data.Curvas_Resf[PositionZ][PosicaoCurvasResfriamento])/
                          (Data.PosInicial[PositionZ+1]- Data.PosInicial[PositionZ]) ) +
                    0.5*( (Data.Curvas_Resf[PositionZ+1][PosicaoCurvasResfriamento+1] -
                           Data.Curvas_Resf[PositionZ][PosicaoCurvasResfriamento+1])/
                          (Data.PosInicial[PositionZ+1]- Data.PosInicial[PositionZ]) );

                dT_dt =
                    0.5*( (Data.Curvas_Resf[PositionZ][PosicaoCurvasResfriamento+1] -
                           Data.Curvas_Resf[PositionZ][PosicaoCurvasResfriamento])/
                          (Data.TimeCurves[PosicaoCurvasResfriamento+1]-
                           Data.TimeCurves[PosicaoCurvasResfriamento]) ) +
                    0.5*( (Data.Curvas_Resf[PositionZ+1][PosicaoCurvasResfriamento+1] -
                           Data.Curvas_Resf[PositionZ+1][PosicaoCurvasResfriamento])/
                          (Data.TimeCurves[PosicaoCurvasResfriamento+1]-
                           Data.TimeCurves[PosicaoCurvasResfriamento]) );


                TTemp = Data.Curvas_Resf[PositionZ][PosicaoCurvasResfriamento]+
                           dT_dz*(yLinha - Data.PosInicial[PositionZ]) +
                           dT_dt*( ti -
                                    Data.TimeCurves[PosicaoCurvasResfriamento]);
                VfP.T = TTemp;
        }


}

double getTThermoCouple(double ti, double yi, struct ThermoCoupleData Data){

        double dT_dz;
        double dT_dt;
        double TTemp;
        int zLinha;
        int PositionZ;
        int PosicaoCurvasResfriamento = 0;
        while(ti > Data.TimeCurves[PosicaoCurvasResfriamento+1])
           PosicaoCurvasResfriamento++;

        if (Data.TimeCurves[PosicaoCurvasResfriamento+1] < 0)
         {
          printf ("Tentando ler mais medidas de T do que o dispon�vel");
         }


         PositionZ = 0;
         while (yi > Data.PosInicial[PositionZ+1])
                PositionZ++;

                dT_dz =
                    0.5*( (Data.Curvas_Resf[PositionZ+1][PosicaoCurvasResfriamento] -
                           Data.Curvas_Resf[PositionZ][PosicaoCurvasResfriamento])/
                          (Data.PosInicial[PositionZ+1]- Data.PosInicial[PositionZ]) ) +
                    0.5*( (Data.Curvas_Resf[PositionZ+1][PosicaoCurvasResfriamento+1] -
                           Data.Curvas_Resf[PositionZ][PosicaoCurvasResfriamento+1])/
                          (Data.PosInicial[PositionZ+1]- Data.PosInicial[PositionZ]) );

                dT_dt =
                    0.5*( (Data.Curvas_Resf[PositionZ][PosicaoCurvasResfriamento+1] -
                           Data.Curvas_Resf[PositionZ][PosicaoCurvasResfriamento])/
                          (Data.TimeCurves[PosicaoCurvasResfriamento+1]-
                           Data.TimeCurves[PosicaoCurvasResfriamento]) ) +
                    0.5*( (Data.Curvas_Resf[PositionZ+1][PosicaoCurvasResfriamento+1] -
                           Data.Curvas_Resf[PositionZ+1][PosicaoCurvasResfriamento])/
                          (Data.TimeCurves[PosicaoCurvasResfriamento+1]-
                           Data.TimeCurves[PosicaoCurvasResfriamento]) );

                TTemp = Data.Curvas_Resf[PositionZ][PosicaoCurvasResfriamento]+
                           dT_dz*(yi - Data.PosInicial[PositionZ]) +
                           dT_dt*( ti -
                                    Data.TimeCurves[PosicaoCurvasResfriamento]);
                return TTemp;

}
#endif

double ThermoCouplePositionT(double y){
        y = Ly - y;
        int iVF;
        double yr;
        double Tinterpol, Tcenter, Tup, Tdown;
        struct MacroMatrix *VFcenter, *VFup,*VFdown;

        iVF = mmay - y/dyVF;
        VFcenter = &VF[iVF];
        if(iVF != 0)
                VFup = &VF[iVF - 1];
        if(iVF != mmacroy - 1)
                VFdown = &VF[iVF + 1];

        yr = y - dyVF*(mmay - iVF - 0.5);

        Tcenter = VFcenter->T;
        Tinterpol = Tcenter;
        if(yr >= 0.0 && iVF != 0){
                Tup = VFup->T;
                Tinterpol += (Tup - Tcenter)*yr/dyVF;
        }
        if(yr <= 0.0 && iVF != mmacroy - 1){
                Tdown = VFdown->T;
                Tinterpol +=(Tcenter - Tdown)*yr/dyVF;
        }
        if(iVF == 0 && yr >= 0.0){
                Tup = Tboundary(Tcenter, Tfarfieldup, (VFcenter->fS)*Ksolid + (1.0 - (VFcenter->fS))*Kliquid, Hconvectionup, dyVF);
                Tinterpol += (Tup - Tcenter)*yr/dyVF;
        }
        if(iVF == mmacroy - 1 && yr <= 0.0){
                Tdown = Tboundary(Tcenter, Tfarfielddown, VFcenter->fS*Ksolid + (1.0 - VFcenter->fS)*Kliquid, Hconvectiondown, dyVF);
                Tinterpol += (Tcenter - Tdown)*yr/dyVF;

        }

        return Tinterpol;

}

double TInterpoly(int interface){

        double Tup, Tdown, Tleft, Tright, Tinterpol;

        double Tcenter;

        int ir, jr; // �ndices do macrosc�pico (a), e �ndices relativos (r)

        double x, y;

        ir = ii%mmiy;
        jr = ji%mmix;

        if( (ii - ir)/mmiy != ia){
                ia=ia;
        }

        if( (ji - jr)/mmix != ja){
                ia=ia;
        }

        /*
        ia = (ii - ir)/mmiy;
        ja = (ji - jr)/mmix;
        */

        Tcenter =  VF[ia*mmax + ja].T;

        Tinterpol = Tcenter;

        x = -dxVF/2.0 + dxCA*(jr + 0.5);
        y = -(-dyVF/2.0 + dyCA*(ir + 0.5));


        if(interface == 0) // SOUTH
                y -= dyCA/2.0;
        if(interface == 1) // NORTH
                y += dyCA/2.0;

        TVF = VfP.T;


        /* Verifica��o se est� na metade da esquerda ou na metade da direita no interior do volume finito */

        if( mmix != 1 && mmax != 1){
        if( jr > mmix/2.0 && ja != mmax - 1){
                Tright = VF[ia*mmax + ja + 1].T;
                Tinterpol += (Tright - Tcenter)*x/dxVF;
        }
        else if( jr <= mmix/2.0 && ja != 0){
                Tleft = VF[ia*mmax + ja - 1].T;
                Tinterpol += (Tcenter - Tleft)*x/dxVF;
        }
        else if( jr <= mmix/2.0 && ja == 0){
                Tleft = Tboundary(Tcenter, Tfarfieldleft, VF[ia*mmax + ja].fS*Ksolid + (1.0 - VF[ia*mmax + ja].fS)*Kliquid, Hconvectionleft, dxVF);
                Tinterpol += (Tcenter - Tleft)*x/dxVF;
        }
        else if( jr > mmiy/2.0 && ja == mmay - 1){
                Tright = Tboundary(Tcenter, Tfarfieldright, VF[ia*mmax + ja].fS*Ksolid + (1.0 - VF[ia*mmax + ja].fS)*Kliquid, Hconvectionright, dxVF);
                Tinterpol += (Tright - Tcenter)*x/dxVF;
        }
        }


        // Verifica��o se est� na metade superior ou na metade inferior do volume finito

        if( mmiy != 1 && mmay != 1){
        if( ir > mmiy/2.0 && ia != mmay - 1){
                Tdown = VF[(ia + 1)*mmax + ja].T;
                Tinterpol +=(Tcenter - Tdown)*y/dyVF;
        }
        else if( ir <= mmiy/2.0 && ia != 0){
                Tup = VF[(ia - 1)*mmax + ja].T;
                Tinterpol += (Tup - Tcenter)*y/dyVF;
        }
        else if( ir <= mmiy/2.0 && ia == 0){
                Tup = Tboundary(Tcenter, Tfarfieldup, VF[ia*mmax + ja].fS*Ksolid + (1.0 - VF[ia*mmax + ja].fS)*Kliquid, Hconvectionup, dyVF);
                Tinterpol += (Tup - Tcenter)*y/dyVF;
        }
        else if( ir > mmiy/2.0 && ia == mmay - 1){
                Tdown = Tboundary(Tcenter, Tfarfielddown, VF[ia*mmax + ja].fS*Ksolid + (1.0 - VF[ia*mmax + ja].fS)*Kliquid, Hconvectiondown, dyVF);
                Tinterpol += (Tcenter - Tdown)*y/dyVF;
        }
        }

        if(Tinterpol < 0.0){
                Tinterpol = Tinterpol;
        }

        if( mmax == 1 && mmay == 1){
                Tinterpol = VF[ia*mmax + ja].T;
        }

        #ifdef NoInterpol
                Tinterpol = VF[ia*mmax + ja].T;
        #endif


        return Tinterpol;
}

double TInterpol(){

        double Tup, Tdown, Tleft, Tright, Tinterpol;

        double Tcenter;

        int ir, jr; // �ndices do macrosc�pico (a), e �ndices relativos (r)

        double x, y;

        ir = ii%mmiy;
        jr = ji%mmix;

        if( (ii - ir)/mmiy != ia){
                ia=ia;
        }

        if( (ji - jr)/mmix != ja){
                ia=ia;
        }

        /*
        ia = (ii - ir)/mmiy;
        ja = (ji - jr)/mmix;
        */

        Tcenter =  VF[ia*mmax + ja].T;

        Tinterpol = Tcenter;

        x = -dxVF/2.0 + dxCA*(jr + 0.5);
        y = -(-dyVF/2.0 + dyCA*(ir + 0.5));

        TVF = VfP.T;


        /* Verifica��o se est� na metade da esquerda ou na metade da direita no interior do volume finito */

        if( mmix != 1 && mmax != 1){
        if( jr > mmix/2.0 && ja != mmax - 1){
                Tright = VF[ia*mmax + ja + 1].T;
                Tinterpol += (Tright - Tcenter)*x/dxVF;
        }
        else if( jr <= mmix/2.0 && ja != 0){
                Tleft = VF[ia*mmax + ja - 1].T;
                Tinterpol += (Tcenter - Tleft)*x/dxVF;
        }
        else if( jr <= mmix/2.0 && ja == 0){
                Tleft = Tboundary(Tcenter, Tfarfieldleft, VF[ia*mmax + ja].fS*Ksolid + (1.0 - VF[ia*mmax + ja].fS)*Kliquid, Hconvectionleft, dxVF);
                Tinterpol += (Tcenter - Tleft)*x/dxVF;
        }
        else if( jr > mmiy/2.0 && ja == mmay - 1){
                Tright = Tboundary(Tcenter, Tfarfieldright, VF[ia*mmax + ja].fS*Ksolid + (1.0 - VF[ia*mmax + ja].fS)*Kliquid, Hconvectionright, dxVF);
                Tinterpol += (Tright - Tcenter)*x/dxVF;
        }
        }


        // Verifica��o se est� na metade superior ou na metade inferior do volume finito

        if( mmiy != 1 && mmay != 1){
        if( ir > mmiy/2.0 && ia != mmay - 1){
                Tdown = VF[(ia + 1)*mmax + ja].T;
                Tinterpol +=(Tcenter - Tdown)*y/dyVF;
        }
        else if( ir <= mmiy/2.0 && ia != 0){
                Tup = VF[(ia - 1)*mmax + ja].T;
                Tinterpol += (Tup - Tcenter)*y/dyVF;
        }
        else if( ir <= mmiy/2.0 && ia == 0){
                Tup = Tboundary(Tcenter, Tfarfieldup, VF[ia*mmax + ja].fS*Ksolid + (1.0 - VF[ia*mmax + ja].fS)*Kliquid, Hconvectionup, dyVF);
                Tinterpol += (Tup - Tcenter)*y/dyVF;
        }
        else if( ir > mmiy/2.0 && ia == mmay - 1){
                Tdown = Tboundary(Tcenter, Tfarfielddown, VF[ia*mmax + ja].fS*Ksolid + (1.0 - VF[ia*mmax + ja].fS)*Kliquid, Hconvectiondown, dyVF);
                Tinterpol += (Tcenter - Tdown)*y/dyVF;
        }
        }

        if(Tinterpol < 0.0){
                Tinterpol = Tinterpol;
        }

        if( mmax == 1 && mmay == 1){
                Tinterpol = VF[ia*mmax + ja].T;
        }

        #ifdef NoInterpol
                Tinterpol = VF[ia*mmax + ja].T;
        #endif


        return Tinterpol;
}


double dHInterpol(){

        double dHup, dHdown, dHleft, dHright;

        double dHcenter, dHinterpol;

        int ir, jr; // �ndices do macrosc�pico, e �ndices relativos

        double x, y;

        ir = ii%mmiy;
        jr = ji%mmix;

        /*
        ia = (ii - ir)/mmiy;
        ja = (ji - jr)/mmix;
        */

        dHcenter =  VF[ia*mmax + ja].dH;

        dHinterpol = dHcenter;

        x = -dxVF/2.0 + dxCA*(jr + 0.5);
        y = -(-dyVF/2.0 + dyCA*(ir + 0.5));



        // Verifica��o se est� na metade da esquerda ou na metade da direita no interior do volume finito

        if( mmax != 1 && mmix != 1){
        if( jr > mmix/2.0 && ja != mmax - 1){
                dHright = VF[ia*mmax + ja + 1].dH;
                dHinterpol += (dHright - dHcenter)*x/dxVF;
        }
        else if( jr <= mmix/2.0 && ja != 0){
                dHleft = VF[ia*mmax + ja - 1].dH;
                dHinterpol += (dHcenter - dHleft)*x/dxVF;
        }
        else if( jr > mmix/2.0 && ja == mmax - 1){
                dHleft = VF[ia*mmax + ja - 1].dH;
                dHinterpol += (dHcenter - dHleft)*x/dxVF;
        }
        else if( jr <= mmix/2.0 && ja == 0){
                dHright = VF[ia*mmax + ja + 1].dH;
                dHinterpol += (dHright - dHcenter)*x/dxVF;
        }
        }


        // Verifica��o se est� na metade superior ou na metade inferior do volume finito

        if( mmay != 1 && mmiy != 1){
        if( ir > mmiy/2.0 && ia != mmay - 1){
                dHdown = VF[(ia + 1)*mmax + ja].dH;
                dHinterpol += (dHcenter - dHdown)*y/dyVF;
        }
        else if( ir <= mmiy/2.0 && ia != 0){
                dHup = VF[(ia - 1)*mmax + ja].dH;
                dHinterpol += (dHup - dHcenter)*y/dyVF;
        }
        else if( ir > mmiy/2.0 && ia == mmay - 1){
                dHup = VF[(ia - 1)*mmax + ja].dH;
                dHinterpol += (dHup - dHcenter)*y/dyVF;
        }
        else if( ir <= mmiy/2.0 && ia == 0){
                dHdown = VF[(ia + 1)*mmax + ja].dH;
                dHinterpol += (dHcenter - dHdown)*y/dyVF;
        }
        }


        if( mmax == 1 && mmay == 1){
                dHinterpol = VF[ia*mmax + ja].dH;
        }

        #ifdef NoInterpol
                dHinterpol = VF[ia*mmax + ja].dH;
        #endif


        return dHinterpol;
}



double fLInterpol(){

        double fLup, fLdown, fLleft, fLright;

        double fLcenter, fLinterpol;

        int ir, jr; // �ndices do macrosc�pico, e �ndices relativos

        double x, y;

        ir = ii%mmiy;
        jr = ji%mmix;

        /*
        ia = (ii - ir)/mmiy;
        ja = (ji - jr)/mmix;
        */

        fLcenter =  VF[ia*mmax + ja].fLCA;

        fLinterpol = fLcenter;

        x = -dxVF/2.0 + dxCA*(jr + 0.5);
        y = -(-dyVF/2.0 + dyCA*(ir + 0.5));



        // Verifica��o se est� na metade da esquerda ou na metade da direita no interior do volume finito

        if( mmax != 1 && mmix != 1){
        if( jr > mmix/2.0 && ja != mmax - 1){
                fLright = VF[ia*mmax + ja + 1].fLCA;
                fLinterpol += (fLright - fLcenter)*x/dxVF;
        }
        else if( jr <= mmix/2.0 && ja != 0){
                fLleft = VF[ia*mmax + ja - 1].fLCA;
                fLinterpol += (fLcenter - fLleft)*x/dxVF;
        }
        else if( jr > mmix/2.0 && ja == mmax - 1){
                fLleft = VF[ia*mmax + ja - 1].fLCA;
                fLinterpol += (fLcenter - fLleft)*x/dxVF;
        }
        else if( jr <= mmix/2.0 && ja == 0){
                fLright = VF[ia*mmax + ja + 1].fLCA;
                fLinterpol += (fLright - fLcenter)*x/dxVF;
        }
        }


        // Verifica��o se est� na metade superior ou na metade inferior do volume finito

        if( mmay != 1 && mmiy != 1){
        if( ir > mmiy/2.0 && ia != mmay - 1){
                fLdown = VF[(ia + 1)*mmax + ja].fLCA;
                fLinterpol += (fLcenter - fLdown)*y/dyVF;
        }
        else if( ir <= mmiy/2.0 && ia != 0){
                fLup = VF[(ia - 1)*mmax + ja].fLCA;
                fLinterpol += (fLup - fLcenter)*y/dyVF;
        }
        else if( ir > mmiy/2.0 && ia == mmay - 1){

        }
        else if( ir <= mmiy/2.0 && ia == 0){

        }
        }


        if( mmax == 1 && mmay == 1){
                fLinterpol = VF[ia*mmax + ja].fLCA;
        }

        #ifdef NoInterpol
                fLinterpol = VF[ia*mmax + ja].fLCA;
        #endif


        return fLinterpol;
}

double edVdVsInterpol(){



        double fLup, fLdown, fLleft, fLright;

        double fLcenter, fLinterpol;

        int ir, jr; // �ndices do macrosc�pico, e �ndices relativos

        double x, y;

        ir = ii%mmiy;
        jr = ji%mmix;

        /*
        ia = (ii - ir)/mmiy;
        ja = (ji - jr)/mmix;
        */

        fLcenter =  (1.0 - VfP.fS - VfP.fE - VfP.fL)*(VfP.Vd - VfP.Vs);

        fLinterpol = fLcenter;

        x = -dxVF/2.0 + dxCA*(jr + 0.5);
        y = -(-dyVF/2.0 + dyCA*(ir + 0.5));



        // Verifica��o se est� na metade da esquerda ou na metade da direita no interior do volume finito

        if( mmax != 1 && mmix != 1){
        if( jr > mmix/2.0 && ja != mmax - 1){
                fLright = (1.0 - VfE.fS - VfE.fE - VfE.fL)*(VfE.Vd - VfE.Vs);
                fLinterpol += (fLright - fLcenter)*x/dxVF;
        }
        else if( jr <= mmix/2.0 && ja != 0){
                fLleft = (1.0 - VfW.fS - VfW.fE - VfW.fL)*(VfW.Vd - VfW.Vs);
                fLinterpol += (fLcenter - fLleft)*x/dxVF;
        }
        else if( jr > mmix/2.0 && ja == mmax - 1){
                fLleft = (1.0 - VfW.fS - VfW.fE - VfW.fL)*(VfW.Vd - VfW.Vs);
                fLinterpol += (fLcenter - fLleft)*x/dxVF;
        }
        else if( jr <= mmix/2.0 && ja == 0){
                fLright = (1.0 - VfE.fS - VfE.fE - VfE.fL)*(VfE.Vd - VfE.Vs);
                fLinterpol += (fLright - fLcenter)*x/dxVF;
        }
        }


        // Verifica��o se est� na metade superior ou na metade inferior do volume finito

        if( mmay != 1 && mmiy != 1){
        if( ir > mmiy/2.0 && ia != mmay - 1){
                fLdown = (1.0 - VfS.fS - VfS.fE - VfS.fL)*(VfS.Vd - VfS.Vs);
                fLinterpol += (fLcenter - fLdown)*y/dyVF;
        }
        else if( ir <= mmiy/2.0 && ia != 0){
                fLup = (1.0 - VfN.fS - VfN.fE - VfN.fL)*(VfN.Vd - VfN.Vs);
                fLinterpol += (fLup - fLcenter)*y/dyVF;
        }
        else if( ir > mmiy/2.0 && ia == mmay - 1){
                
        }
        else if( ir <= mmiy/2.0 && ia == 0){
                
        }
        }


        if( mmax == 1 && mmay == 1){
                fLinterpol = (1.0 - VfP.fS - VfP.fE - VfP.fL)*(VfP.Vd - VfP.Vs);
        }

        #ifdef NoInterpol
                fLinterpol = (1.0 - VfP.fS - VfP.fE - VfP.fL)*(VfP.Vd - VfP.Vs);
        #endif


        return fLinterpol;
}




double NeInterpol(){



        double fLup, fLdown, fLleft, fLright;

        double fLcenter, fLinterpol;

        int ir, jr; // �ndices do macrosc�pico, e �ndices relativos

        double x, y;

        ir = ii%mmiy;
        jr = ji%mmix;

        /*
        ia = (ii - ir)/mmiy;
        ja = (ji - jr)/mmix;
        */

        fLcenter =  VfP.Ne;

        fLinterpol = fLcenter;

        x = -dxVF/2.0 + dxCA*(jr + 0.5);
        y = -(-dyVF/2.0 + dyCA*(ir + 0.5));



        // Verifica��o se est� na metade da esquerda ou na metade da direita no interior do volume finito

        if( mmax != 1 && mmix != 1){
        if( jr > mmix/2.0 && ja != mmax - 1){
                fLright = VfE.Ne;
                fLinterpol += (fLright - fLcenter)*x/dxVF;
        }
        else if( jr <= mmix/2.0 && ja != 0){
                fLleft = VfW.Ne;
                fLinterpol += (fLcenter - fLleft)*x/dxVF;
        }
        else if( jr > mmix/2.0 && ja == mmax - 1){
                fLleft = VfW.Ne;
                fLinterpol += (fLcenter - fLleft)*x/dxVF;
        }
        else if( jr <= mmix/2.0 && ja == 0){
                fLright = VfE.Ne;
                fLinterpol += (fLright - fLcenter)*x/dxVF;
        }
        }


        // Verifica��o se est� na metade superior ou na metade inferior do volume finito

        if( mmay != 1 && mmiy != 1){
        if( ir > mmiy/2.0 && ia != mmay - 1){
                fLdown = VfS.Ne;
                fLinterpol += (fLcenter - fLdown)*y/dyVF;
        }
        else if( ir <= mmiy/2.0 && ia != 0){
                fLup = VfN.Ne;
                fLinterpol += (fLup - fLcenter)*y/dyVF;
        }
        else if( ir > mmiy/2.0 && ia == mmay - 1){

        }
        else if( ir <= mmiy/2.0 && ia == 0){

        }
        }


        if( mmax == 1 && mmay == 1){
                fLinterpol = VfP.Ne;
        }

        #ifdef NoInterpol
                fLinterpol = VfP.Ne;
        #endif


        return fLinterpol;
}




double CLInterpol(){

        double CLup, CLdown, CLleft, CLright;

        double CLcenter, CLinterpol;

        int ir, jr; // �ndices do macrosc�pico, e �ndices relativos

        double x, y;

        ir = ii%mmiy;
        jr = ji%mmix;

        /*
        ia = (ii - ir)/mmiy;
        ja = (ji - jr)/mmix;
        */

        CLcenter =  VF[ia*mmax + ja].CL;

        CLinterpol = CLcenter;

        x = -dxVF/2.0 + dxCA*(jr + 0.5);
        y = -(-dyVF/2.0 + dyCA*(ir + 0.5));



        // Verifica��o se est� na metade da esquerda ou na metade da direita no interior do volume finito

        if( mmax != 1 && mmix != 1){
        if( jr > mmix/2.0 && ja != mmax - 1){
                CLright = VF[ia*mmax + ja + 1].CL;
                CLinterpol += (CLright - CLcenter)*x/dxVF;
        }
        else if( jr <= mmix/2.0 && ja != 0){
                CLleft = VF[ia*mmax + ja - 1].CL;
                CLinterpol += (CLcenter - CLleft)*x/dxVF;
        }
        else if( jr > mmix/2.0 && ja == mmax - 1){
                CLleft = VF[ia*mmax + ja - 1].CL;
                CLinterpol += (CLcenter - CLleft)*x/dxVF;
        }
        else if( jr <= mmix/2.0 && ja == 0){
                CLright = VF[ia*mmax + ja + 1].CL;
                CLinterpol += (CLright - CLcenter)*x/dxVF;
        }
        }


        // Verifica��o se est� na metade superior ou na metade inferior do volume finito

        if( mmay != 1 && mmiy != 1){
        if( ir > mmiy/2.0 && ia != mmay - 1){
                CLdown = VF[(ia + 1)*mmax + ja].CL;
                CLinterpol += (CLcenter - CLdown)*y/dyVF;
        }
        else if( ir <= mmiy/2.0 && ia != 0){
                CLup = VF[(ia - 1)*mmax + ja].CL;
                CLinterpol += (CLup - CLcenter)*y/dyVF;
        }
        else if( ir > mmiy/2.0 && ia == mmay - 1){
                CLup = VF[(ia - 1)*mmax + ja].CL;
                CLinterpol += (CLup - CLcenter)*y/dyVF;
        }
        else if( ir <= mmiy/2.0 && ia == 0){
                CLdown = VF[(ia + 1)*mmax + ja].CL;
                CLinterpol += (CLcenter - CLdown)*y/dyVF;
        }
        }


        if( mmax == 1 && mmay == 1){
                CLinterpol = VF[ia*mmax + ja].CL;
        }




        return CLinterpol;
}




void TipInterpol(){

        


}

