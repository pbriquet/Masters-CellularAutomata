#include <math.h>
#include <stdio.h>
#include <stdlib.h>

//---------------------------------------------------------------------------
#pragma hdrstop

#include "GrowthKnetics.h"
//---------------------------------------------------------------------------
#pragma package(smart_init)



#include "GlobalDefinitions.h"
#include "SetupDatabase.h"
#include "GlobalVariables.h"


#include "ExtraFunctions.h"
#include "InterpolationFunctions.h"

#ifdef eglocal
void CalculateAmax(){

                double dr[4][2];           // Distância entre o centro do envelope e o centro da célula a ser ativada
                double n[4][2];         // Vetor unitário que direciona os semi-lados do envelope. O primeiro índice: 0 = vetor 11, 1 = vetor I1, 2 = vetor II, 3 = vetor 1I. Segundo índice: 0 = x, 1 = y
                double costheta;
                double sintheta;

                double AreaMax = 0.0;
                double Area;
                double sumN[2];
                double Lmax = 0.0;

                for( int k=0 ; k<=3; k++){

                dr[k][0]= ( (k-1)%2 )*(-dxCA) - CA[ii*mx + ji].dC[0];       // Obtém dr[k] do sítio vizinho selecionado
                dr[k][1]= ( (k-2)%2 )*(-dyCA) - CA[ii*mx + ji].dC[1];
                }

                costheta = cos(CA[ii*mx + ji].theta );
                sintheta = sin(CA[ii*mx + ji].theta );

                n[0][0] = costheta;
                n[0][1] = sintheta;
                n[1][0] = -sintheta;
                n[1][1] = costheta;
                n[2][0] = -costheta;
                n[2][1] = -sintheta;
                n[3][0] = sintheta;
                n[3][1] = -costheta;

                for( int k=0; k<=3 ; k++){

                                for(int m=0; m<=3; m++){
                                                sumN[0] = n[m][0];
                                                if(m!=3)
                                                sumN[0] += n[m+1][0];
                                                else
                                                sumN[0] += n[0][0];

                                                sumN[1] = n[m][1];
                                                if(m!=3)
                                                sumN[1] += n[m+1][1];
                                                else
                                                sumN[1] += n[0][1];

                                                if( proj(dr[k],sumN) >= 0.0){
                                                Lmax = proj(dr[k], sumN);
                                                Area = 2.0*Lmax*Lmax;

                                                if(Area > 0.0 && Area > AreaMax)
                                                                AreaMax = Area;
                                                }
                                }

                }

                CA[ii*mx + ji].Amax = AreaMax;



                

}
#endif

double ThermalVelocity(double dT){
        if(dT > 0){
                #ifndef THERMALVELOCITY2
                return a2*pow(dT,2) + a3*pow(dT,3);
                //return 1.0e-4;
                #endif
                #ifdef THERMALVELOCITY2
                return 3.0e-6*pow(dT , 2.7);
                #endif
        }
        else
                return 0;
}

double SolutalVelocity( double Cd, double Cl){

                #ifdef DendriticTipEquilibriumWithVFTemperature
                                //Cd = (VF[ia*mmax + ja].T - Tf)/ml;
                                Cd = (TInterpol() - Tf)/ml;
                #endif
        double omega = (Cd - Cl)/Cd/(1.0-kpart);
        Cd = Cd;
        if( Cl >= Cd )
                return 0.0;

        return DiffLiq*ml*(kpart - 1.0)*Cd/M_PI/M_PI/GibbsThompson*pow(0.4567*pow(omega/(1.0 - omega) , 1.195), 2.0);

}

double DendriticVelocity(){

        double V = 0.0;

        #ifdef Cofarfield
                CLVFTipGrowth = Co;
        #endif

        #ifdef CLVFTipGrowthfarfield
                CLVFTipGrowth = VF[ia*mmax + ja].CL;
        #endif

        #ifdef CLInterpolfarfield
                CLVFTipGrowth = CLInterpol();
        #endif

        #ifdef SOLUTALVELOCITY
                V = SolutalVelocity( CaP.Cd, CLVFTipGrowth);
        #endif

        #ifdef THERMALVELOCITY
                V = ThermalVelocity( Tliq - TCA);
        #endif

        #ifdef CONSTANTVELOCITY
                V = 0.0001;
        #endif



        return V;

}






void CAGrowth(){





        int l;
        ii=ii;

        // 3.4) Variáveis do algoritmo de crescimento


        int criteria = 0;        // Contador de projeções que satisfazem o critério, o que verifica se a célula foi ativada, cada vetor é para contar lados opostos


        double dL[4];              // Variação de L do sítio


        int desactivator;       // Contador de sítios sólidos em volta de uma célula



                // 3.4.1) Vetores para Crescimento e Ativação


        double dr[4][2];           // Distância entre o centro do envelope e o centro da célula a ser ativada
        double n[4][2];         // Vetor unitário que direciona os semi-lados do envelope. O primeiro índice: 0 = vetor 11, 1 = vetor I1, 2 = vetor II, 3 = vetor 1I. Segundo índice: 0 = x, 1 = y

        double costheta;
        double sintheta;

        double modS1;
        double modS2;
        double alfa1;   // Parâmetro alfa para verificação de captura
        double alfa2;
        double beta;    // Parâmetro para cálculo das diagonais do novo polígono de crescimento.
        double R;

        // Armazenamentos temporários para cálculo na captura.
        double modulodif;       // O módulo da subtração dos vetores diagonais
        double prodescalar1;    // Produto escalar de dr com a primeira diagonal
        double prodescalar2;    // Produto escalar de dr com a segunda diagonal
        double IS1, IS2;





        double who[4];          // Memoriza quem ativou: [0][0] - Direita, [1][0] - De cima, [0][1] - Esquerda, [1][1] - Embaixo, e a projeção que atendeu ao critério

        for(l=0;l<=3;l++){
                who[l]=-1;      // "Zera" o "who"
        }


        double Lnovo;

        double dCn1[2];          // Memoriza descentralização para o novo envelope
        double dCn2[2];


                // 3.4.2) Ponteiros para vizinhos da célula

        double *neighborstheta[4];
        double *neighborsL[4][4];
        double *neighborsdC[2][4];
        int *neighborsstate[4];
        int *neighborCd[4];

        struct MicroMatrix *neighbor[4];

        for(int k = 0; k < 4; k++){

        neighborstheta[k] = NULL;
        neighborsL[0][k] = NULL;
        neighborsL[1][k] = NULL;
        neighborsL[2][k] = NULL;
        neighborsL[3][k] = NULL;
        neighborsdC[0][k] = NULL;
        neighborsdC[1][k] = NULL;
        neighborsstate[k] = NULL;
        neighborCd[k]= NULL;

        neighbor[k] = NULL;

        }



                if( Tliq - TCA >= 0.0){

                /*
                dL[0] = getV( Tliquidus - TCA, modulo( clusters[ clusters[ (int) fabs( CA[ii*mx + ji].state ) - 2].cluster].V ), gettheta( clusters[ clusters[ (int) fabs( CA[ii*mx + ji].state ) - 2].cluster].V, CA[ii*mx + ji].theta, 0, t), Co)*dt;
                dL[1] = getV( Tliquidus - TCA, modulo( clusters[ clusters[ (int) fabs( CA[ii*mx + ji].state ) - 2].cluster].V ), gettheta( clusters[ clusters[ (int) fabs( CA[ii*mx + ji].state ) - 2].cluster].V, CA[ii*mx + ji].theta, 1, t), Co)*dt;
                dL[2] = getV( Tliquidus - TCA, modulo( clusters[ clusters[ (int) fabs( CA[ii*mx + ji].state ) - 2].cluster].V ), gettheta( clusters[ clusters[ (int) fabs( CA[ii*mx + ji].state ) - 2].cluster].V, CA[ii*mx + ji].theta, 2, t), Co)*dt;
                dL[3] = getV( Tliquidus - TCA, modulo( clusters[ clusters[ (int) fabs( CA[ii*mx + ji].state ) - 2].cluster].V ), gettheta( clusters[ clusters[ (int) fabs( CA[ii*mx + ji].state ) - 2].cluster].V, CA[ii*mx + ji].theta, 3, t), Co)*dt;
                */

                #ifndef NOGrowth
                dL[0] = DendriticVelocity()*dt;
                #endif
                #ifdef NOGrowth
                dL[0] = 0.0;
                #endif
                //dL[0] = SolutalVelocity(CaP.Cd , VfP.CL )*dt;
                //dL[0] = ThermalVelocity( Tliquidus - TCA)*dt;
                dL[1] = dL[0];
                dL[2] = dL[0];
                dL[3] = dL[0];

                if( dL[0] > dyCA )
                        printf("\n\nAlert! CA Growth Criteria.");

                #ifdef eglocal
                fGCA = ( 2.0*CaP.L[0]*CaP.L[0] - CaP.Amin)/(CaP.Amax - CaP.Amin);


                if( fGCA > 1.0 )
                                ii=ii;

                #endif


                /*
                dL[0] = V( Tliquidus - TCA)*dt;
                dL[1] = V( Tliquidus - TCA)*dt;
                dL[2] = V( Tliquidus - TCA)*dt;
                dL[3] = V( Tliquidus - TCA)*dt;
                */

                /*
                dL[0] = 0.0;
                dL[1] = 0.0;
                dL[2] = 0.0;
                dL[3] = 0.0;
                */

                CA[ii*mx + ji].L[0] += dL[0];   // Cresce sólidos ativos multiplicando a velocidade da dendrita na diagonal do envelope com dt
                CA[ii*mx + ji].L[1] += dL[1];
                CA[ii*mx + ji].L[2] += dL[2];
                CA[ii*mx + ji].L[3] += dL[3];



                #ifdef eglocal

                fGCAdt = ( 2.0*CaP.L[0]*CaP.L[0] - CaP.Amin)/(CaP.Amax - CaP.Amin);

                CA[ii*mx + ji].dfG = fGCAdt - fGCA;

                if( CA[ii*mx + ji].dfG > 0.5)
                                ii=ii;

                #endif



                // Verificação de Ativação dos sítios vizinhos


                if(ii != my - 1 && CA[(ii+1)*mx + ji].state==0){  // Verifica se sítio está na extremo inferior da malha, e se não está e sítio inferior é líquido, aponta para o sítio
                neighborsL[0][3] = &CA[(ii+1)*mx + ji].L[0];
                neighborsL[1][3] = &CA[(ii+1)*mx + ji].L[1];
                neighborsL[2][3] = &CA[(ii+1)*mx + ji].L[2];
                neighborsL[3][3] = &CA[(ii+1)*mx + ji].L[3];
                neighborsstate[3] = &CA[(ii+1)*mx + ji].state;       // k = 0 -> Direita ; k = 1 -> Superior ; k = 2 -> Esquerda ; k = 3 -> Inferior
                neighborstheta[3] = &CA[(ii+1)*mx + ji].theta;
                neighborsdC[0][3] = &CA[(ii+1)*mx + ji].dC[0];
                neighborsdC[1][3] = &CA[(ii+1)*mx + ji].dC[1];

                neighbor[3] = &CA[(ii+1)*mx + ji];

                }

                if(ji != mx - 1 &&  CA[ii*mx + ji + 1].state==0){  // Verifica se sítio está na extrema direita da malha, e se não está e sítio direito é líquido, aponta para o sítio
                neighborsL[0][0] = &CA[ii*mx + ji + 1].L[0];
                neighborsL[1][0] = &CA[ii*mx + ji + 1].L[1];
                neighborsL[2][0] = &CA[ii*mx + ji + 1].L[2];
                neighborsL[3][0] = &CA[ii*mx + ji + 1].L[3];
                neighborsstate[0] = &CA[ii*mx + ji + 1].state;
                neighborstheta[0] = &CA[ii*mx + ji + 1].theta;
                neighborsdC[0][0] = &CA[ii*mx + ji + 1].dC[0];
                neighborsdC[1][0] = &CA[ii*mx + ji + 1].dC[1];

                neighbor[0] = &CA[ii*mx + ji + 1];
                }

                if(ii != 0 && CA[(ii-1)*mx + ji].state==0){   // Verifica se sítio está na extremo superior da malha, e se não está e sítio superior é líquido, aponta para o sítio
                neighborsL[0][1] = &CA[(ii-1)*mx + ji].L[0];
                neighborsL[1][1] = &CA[(ii-1)*mx + ji].L[1];
                neighborsL[2][1] = &CA[(ii-1)*mx + ji].L[2];
                neighborsL[3][1] = &CA[(ii-1)*mx + ji].L[3];
                neighborsstate[1] = &CA[(ii-1)*mx + ji].state;
                neighborstheta[1] = &CA[(ii-1)*mx + ji].theta;
                neighborsdC[0][1] = &CA[(ii-1)*mx + ji].dC[0];
                neighborsdC[1][1] = &CA[(ii-1)*mx + ji].dC[1];

                neighbor[1] = &CA[(ii-1)*mx + ji];
                }

                if(ji != 0 && CA[ii*mx + ji - 1].state==0){   // Verifica se sítio está na extrema esquerda da malha, e se não está e sítio esquerdo é líquido, aponta para o sítio
                neighborsL[0][2] = &CA[ii*mx + ji - 1].L[0];
                neighborsL[1][2] = &CA[ii*mx + ji - 1].L[1];
                neighborsL[2][2] = &CA[ii*mx + ji - 1].L[2];
                neighborsL[3][2] = &CA[ii*mx + ji - 1].L[3];
                neighborsstate[2] = &CA[ii*mx + ji - 1].state;
                neighborstheta[2] = &CA[ii*mx + ji - 1].theta;
                neighborsdC[0][2] = &CA[ii*mx + ji - 1].dC[0];
                neighborsdC[1][2] = &CA[ii*mx + ji - 1].dC[1];

                neighbor[2] = &CA[ii*mx + ji - 1];
                }



                                                                // Cálculo dos vetores para verificação de critério

                costheta = cos(CA[ii*mx + ji].theta );
                sintheta = sin(CA[ii*mx + ji].theta );



                n[0][0] = costheta;
                n[0][1] = sintheta;
                n[1][0] = -sintheta;
                n[1][1] = costheta;
                n[2][0] = -costheta;
                n[2][1] = -sintheta;
                n[3][0] = sintheta;
                n[3][1] = -costheta;

                for( int k=0 ; k<=3; k++){

                dr[k][0]= ( (k-1)%2 )*(-dxCA) - CA[ii*mx + ji].dC[0];       // Obtém dr[k] do sítio vizinho selecionado
                dr[k][1]= ( (k-2)%2 )*(-dyCA) - CA[ii*mx + ji].dC[1];
                }




                for(int k=0; k<=3 && CA[ii*mx + ji].state >=2 ;k++){
                if(neighborsstate[k]==NULL){     // Sítio vizinho é líquido ? Se não for -> Continue
                continue;
                }



                // Critério de Ativação              // Posso mudar para verificação a verificação do lado que ativou através da verificação se lado atende ao critério com e não atende sem dL

                if( CA[ii*mx + ji].L[0] != 0.0 && CA[ii*mx + ji].L[1] != 0.0 && criteria == 0 ){

                prodescalar1 = proj(dr[k], n[0]);
                alfa1 = prodescalar1/CA[ii*mx +ji].L[0];
                if( alfa1 >= 0.0 && alfa1 <= 1.0){

                prodescalar2 = proj(dr[k], n[1]);
                alfa2 = prodescalar2/CA[ii*mx +ji].L[1];

                if( alfa2 >= 0.0 && alfa2 < 1.0 - alfa1){
                modulodif = sqrt( pow(CA[ii*mx + ji].L[0],2.0) + pow(CA[ii*mx + ji].L[1],2.0) );

                criteria = 1;
                modS1 = CA[ii*mx + ji].L[0];
                modS2 = CA[ii*mx + ji].L[1];

                }

                }
                }

                if( CA[ii*mx + ji].L[1] != 0.0 && CA[ii*mx + ji].L[2] != 0.0 && criteria == 0 ){

                prodescalar1 = proj(dr[k], n[1]);
                alfa1 = prodescalar1/CA[ii*mx +ji].L[1];
                if( alfa1 >= 0.0 && alfa1 <= 1.0){

                prodescalar2 = proj( dr[k], n[2] );
                alfa2 = prodescalar2/CA[ii*mx +ji].L[2];

                if( alfa2 >= 0.0 && alfa2 < 1.0 - alfa1){
                modulodif = sqrt( pow(CA[ii*mx + ji].L[1],2.0) + pow(CA[ii*mx + ji].L[2],2.0) );

                criteria = 2;
                modS1 = CA[ii*mx + ji].L[1];
                modS2 = CA[ii*mx + ji].L[2];

                }

                }
                }

                if( CA[ii*mx + ji].L[2] != 0.0 && CA[ii*mx + ji].L[3] != 0.0 && criteria == 0 ){

                prodescalar1 = proj(dr[k], n[2]);
                alfa1 = prodescalar1/CA[ii*mx +ji].L[2];
                if( alfa1 >= 0.0 && alfa1 <= 1.0){

                prodescalar2 = proj( dr[k], n[3] );
                alfa2 = prodescalar2/CA[ii*mx +ji].L[3];

                if( alfa2 >= 0.0 && alfa2 < 1.0 - alfa1){
                modulodif = sqrt( pow(CA[ii*mx + ji].L[2],2.0) + pow(CA[ii*mx + ji].L[3],2.0) );

                criteria = 3;
                modS1 = CA[ii*mx + ji].L[2];
                modS2 = CA[ii*mx + ji].L[3];

                }

                }
                }

                if( CA[ii*mx + ji].L[3] != 0.0 && CA[ii*mx + ji].L[0] != 0.0 && criteria == 0 ){

                prodescalar1 = proj(dr[k], n[3]);
                alfa1 = prodescalar1/CA[ii*mx +ji].L[3];
                if( alfa1 >= 0 && alfa1 <= 1){

                prodescalar2 = proj( dr[k], n[0] );
                alfa2 = prodescalar2/CA[ii*mx +ji].L[0];

                if( alfa2 >= 0.0 && alfa2 < 1.0 - alfa1){
                modulodif = sqrt( pow(CA[ii*mx + ji].L[3],2.0) + pow(CA[ii*mx + ji].L[0],2.0) );

                criteria = 4;
                modS1 = CA[ii*mx + ji].L[3];
                modS2 = CA[ii*mx + ji].L[0];

                }

                }
                }




                if(criteria > 0 ){       // Se pelo menos 1 dos lados inferior/superior atende ao critério e (AND) pelo menos 1 dos lados esquerdo/direito atende ao critério, então sítio foi capturado

                if(modulodif <= 0.00001){
                modulodif = modulodif;
                }
                beta = ( prodescalar2*modS2 - prodescalar1*modS1 + pow( modS1 , 2.0) )/pow( modulodif , 2.0);

                IS1 = beta*modulodif;
                IS2 = ( 1 - beta )*modulodif;

                Lnovo = (min( IS1, sqrt(2.0)*dxCA ) + min( IS2, sqrt(2.0)*dxCA) )/2.0;

                R = Lnovo/modulodif;

                // Calculando L novo
                for( int m = 0; m < 4; m++){
                neighbor[k]->L[m] = R*CA[ii*mx + ji].L[m];
                }
                neighbor[k]->theta = CA[ii*mx + ji].theta + 10.0;
                neighbor[k]->state = CA[ii*mx + ji].state;        // Sólido Recém-Ativo

                #ifdef CdlEqualFluxinDiffusion
                fdP = 1.0 - CaP.fS - CaP.fE;
                velocity = DendriticVelocity();
                CLVFDiffusion = VfP.CL;
                CdCA = CaP.Cd;
                #ifdef SeCACorrectionSphere
                velocity = 1.0/sqrt(2.0*MPI);
                #endif
                DelAdjustmentFunction();
                Rdl = ExtraliquidDiffusion/del_dl;
                Rp = InterliquidDiffusion/dxCA*2.0*fdP;
                CdlCA = (Rp*CdCA + Rdl*CLVFDiffusion)/(Rp + Rdl);
                #endif

                neighbor[k]->Cd = Cdl;

                #ifdef eglocal
                //neighbor[k]->Amin = (neighbor[k]->L[0] + neighbor[k]->L[2])*(neighbor[k]->L[1] + neighbor[k]->L[3])/2.0;
                neighbor[k]->Amin = (2.0*neighbor[k]->L[0]*neighbor[k]->L[0]);
                #endif

                                if( criteria == 1){

                                                                                        dCn1[0] = (1 - R)*CA[ii*mx +ji].L[0]*n[0][0] - dr[k][0];
                                                                                        dCn1[1] = (1 - R)*CA[ii*mx +ji].L[0]*n[0][1] - dr[k][1];

                                                                                        dCn2[0] = (1 - R)*CA[ii*mx +ji].L[1]*n[1][0] - dr[k][0];
                                                                                        dCn2[1] = (1 - R)*CA[ii*mx +ji].L[1]*n[1][1] - dr[k][1];

                                                                                        if( pow( dCn1[0], 2.0) + pow( dCn1[1], 2.0) <= pow( dCn2[0] , 2.0) + pow( dCn2[1], 2.0)  ){

                                                                                                *neighborsdC[0][k] = dCn1[0];
                                                                                                *neighborsdC[1][k] = dCn1[1];

                                                                                                //*neighborsL[2][k] = 0.0;

                                                                                        }

                                                                                        else{

                                                                                                *neighborsdC[0][k] = dCn2[0];
                                                                                                *neighborsdC[1][k] = dCn2[1];

                                                                                                //*neighborsL[3][k] = 0.0;
                                                                                        }



                                                                                }

                                                                                if( criteria == 2){

                                                                                        dCn1[0] = (1 - R)*CA[ii*mx +ji].L[1]*n[1][0] - dr[k][0];
                                                                                        dCn1[1] = (1 - R)*CA[ii*mx +ji].L[1]*n[1][1] - dr[k][1];

                                                                                        dCn2[0] = (1 - R)*CA[ii*mx +ji].L[2]*n[2][0] - dr[k][0];
                                                                                        dCn2[1] = (1 - R)*CA[ii*mx +ji].L[2]*n[2][1] - dr[k][1];

                                                                                        if( pow( dCn1[0], 2.0) + pow( dCn1[1], 2.0) <= pow( dCn2[0] , 2.0) + pow( dCn2[1], 2.0)  ){

                                                                                                *neighborsdC[0][k] = dCn1[0];
                                                                                                *neighborsdC[1][k] = dCn1[1];

                                                                                                //*neighborsL[3][k] = 0.0;

                                                                                        }

                                                                                        else{

                                                                                                *neighborsdC[0][k] = dCn2[0];
                                                                                                *neighborsdC[1][k] = dCn2[1];

                                                                                                //*neighborsL[0][k] = 0.0;

                                                                                        }

                                                                                }

                                                                                if( criteria == 3){

                                                                                        dCn1[0] = (1 - R)*CA[ii*mx +ji].L[2]*n[2][0] - dr[k][0];
                                                                                        dCn1[1] = (1 - R)*CA[ii*mx +ji].L[2]*n[2][1] - dr[k][1];

                                                                                        dCn2[0] = (1 - R)*CA[ii*mx +ji].L[3]*n[3][0] - dr[k][0];
                                                                                        dCn2[1] = (1 - R)*CA[ii*mx +ji].L[3]*n[3][1] - dr[k][1];

                                                                                        if( pow( dCn1[0], 2.0) + pow( dCn1[1], 2.0) <= pow( dCn2[0] , 2.0) + pow( dCn2[1], 2.0)  ){

                                                                                                *neighborsdC[0][k] = dCn1[0];
                                                                                                *neighborsdC[1][k] = dCn1[1];

                                                                                                //*neighborsL[0][k] = 0.0;

                                                                                        }

                                                                                        else{

                                                                                                *neighborsdC[0][k] = dCn2[0];
                                                                                                *neighborsdC[1][k] = dCn2[1];

                                                                                                //*neighborsL[1][k] = 0.0;

                                                                                        }

                                                                                }

                                                                                if( criteria == 4){

                                                                                        dCn1[0] = (1 - R)*CA[ii*mx +ji].L[3]*n[3][0] - dr[k][0];
                                                                                        dCn1[1] = (1 - R)*CA[ii*mx +ji].L[3]*n[3][1] - dr[k][1];

                                                                                        dCn2[0] = (1 - R)*CA[ii*mx +ji].L[0]*n[0][0] - dr[k][0];
                                                                                        dCn2[1] = (1 - R)*CA[ii*mx +ji].L[0]*n[0][1] - dr[k][1];

                                                                                        if( pow( dCn1[0], 2.0) + pow( dCn1[1], 2.0) <= pow( dCn2[0] , 2.0) + pow( dCn2[1], 2.0)  ){

                                                                                                *neighborsdC[0][k] = dCn1[0];
                                                                                                *neighborsdC[1][k] = dCn1[1];

                                                                                                //*neighborsL[1][k] = 0.0;

                                                                                        }

                                                                                        else{

                                                                                                *neighborsdC[0][k] = dCn2[0];
                                                                                                *neighborsdC[1][k] = dCn2[1];

                                                                                                //*neighborsL[2][k] = 0.0;

                                                                                        }

                                                                                }



                                                                        } // Fim do if(criteria)

                                                                        criteria = 0; // Zera contador do critério de ativação



                                                                } // Fim do for(k)

                                                                desactivator=0;


                                                        // Verificação para desativar célula em crescimento

                                                                #ifndef eglocal
                                                                if(ii == (my - 1) || CA[(ii+1)*mx + ji].state!=0)
                                                                        desactivator++;
                                                                if(ii == 0 || CA[(ii-1)*mx + ji].state!=0)
                                                                        desactivator++;
                                                                if(ji == (mx - 1) || CA[ii*mx + ji + 1].state!=0)
                                                                        desactivator++;
                                                                if(ji == 0 || CA[ii*mx + ji - 1].state!=0)
                                                                        desactivator++;
                                                                #endif

                                                                #ifdef eglocal
                                                                if( (2.0*CaP.L[0]*CaP.L[0]) >= CaP.Amax){
                                                                                fGCAdt = 1.0;
                                                                                CA[ii*mx + ji].L[0] = sqrt( (fGCAdt*(CaP.Amax - CaP.Amin) + CaP.Amin)/2.0);
                                                                                CA[ii*mx + ji].L[1] = CA[ii*mx + ji].L[0];
                                                                                CA[ii*mx + ji].L[2] = CA[ii*mx + ji].L[0];
                                                                                CA[ii*mx + ji].L[3] = CA[ii*mx + ji].L[0];
                                                                                CA[ii*mx + ji].dfG = fGCAdt - fGCA;
                                                                                if(CA[ii*mx + ji].dfG < 0.0)
                                                                                                ii=ii;
                                                                                desactivator = 4;
                                                                }
                                                                #endif



                                                                if(desactivator==4){
                                                                        CA[ii*mx + ji].state= -CA[ii*mx + ji].state;     // Mostra que é um sólido interno // Desativa o sólido em crescimento
                                                                }




                                                        // Verificação para desativar célula em crescimento

                                                                for( int k=0;k<4;k++){       // "Zera" o ponteiro neighbors
                                                                                neighborsL[0][k] = NULL;
                                                                                neighborsL[1][k] = NULL;
                                                                                neighborsL[2][k] = NULL;
                                                                                neighborsL[3][k] = NULL;
                                                                                neighborsstate[k] = NULL;
                                                                                neighborstheta[k] = NULL;
                                                                                neighborsdC[0][k] = NULL;
                                                                                neighborsdC[1][k] = NULL;

                                                                                neighbor[k] = NULL;
                                                                }


                                                        } // Fim do if( Tliquidus - TCA)



}










