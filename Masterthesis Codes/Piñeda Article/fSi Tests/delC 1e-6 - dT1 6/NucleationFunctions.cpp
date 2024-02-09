
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

//---------------------------------------------------------------------------
#pragma hdrstop

#include "NucleationFunctions.h"
//---------------------------------------------------------------------------
#pragma package(smart_init)


#include "GlobalDefinitions.h"
#include "SetupDatabase.h"
#include "GlobalVariables.h"

#include "InterpolationFunctions.h"


void Nucleation(){


                if( CaP.dT >= 0.0 && (Tliq - TCA) >= CaP.dT && CaP.state == 0){

                                int alea = 0;
                                if(ii == 261 && ji == 13)
                                        ii=ii;
                                CaP.L[0] = Lmin;
                                CaP.L[1] = 0.0;
                                CaP.L[2] = 0.0;
                                CaP.L[3] = 0.0;

                                alea = random(nclassestheta - 1);
                                CaP.theta = -M_PI/4.0 + alea/((double) nclassestheta - 1.0)*M_PI/2.0 ;       // Sorteia uma classe dentre as 48 entre [-45o, 45o]
                                NTotalGrains = NTotalGrains + 1;                                                                            // O + 10.0 é para não haver transformação de fase ao nuclear

                                CaP.state = NTotalGrains + 1;
                                CaP.dT = -1.0;

                                CdlCA = VfP.CL;

                                CaP.Cd = Cdl;

                                #ifdef eglocal
                                CaP.Amin = 0.0;
                                #endif

                                clusters[NTotalGrains - 1].thetaclass = alea;
                                clusters[NTotalGrains - 1].cluster = NTotalGrains - 1;
                                clusters[NTotalGrains - 1].size = 1;

                                #ifdef esiCluster
                                        clusters[NTotalGrains - 1].fS = fSmax;
                                        clusters[NTotalGrains - 1].de = sqrt( 2.0 )*sqrt(dyCA*dxCA)/pow( M_PI, 1.0/3.0 );
                                        clusters[NTotalGrains - 1].Ve = M_PI*pow(clusters[NTotalGrains - 1].de, 3.0)/6.0;
                                #endif

                                CaP.theta -= 10.0;


                }



}



























double normalnucS(double dT){   // dn/d(dT) para nucleação na superfície
        return  nsmax*exp(-(pow( (dT - dTsnuc)/dTssigma,2.0)/2.0) )/(sqrt(2.0*M_PI)*dTssigma);
}

                // 2.9.2) dn/d(dT) para nucleação no volume V

double normalnucV(double dT){   // dn/d(dT) para nucleação no volume
        #ifdef NormalDistribution
        return  nvmax*exp(-(pow( (dT - dTvnuc)/dTvsigma,2.0)/2.0) )/(sqrt(2.0*M_PI)*dTvsigma);
        #endif
        #ifdef LogNormalDistribution
        if( dT == 0.0 )
                return 0.0;
        else
                return nvmax/ShapeParameter/sqrt(2.0*M_PI)/dT*exp(-1.0/2.0*pow( log(dT/dTvnuc)/ShapeParameter, 2.0) );
        #endif
}

                // 2.9.3) Fatia de área obtida da função normalnucS pelo Método dos Trapézios

double trapezionormalnucS(double dT){   // Integral pelo método dos trapézios, de dn/d(dT) para nucleação na superfície
        return  (normalnucS(dT) + normalnucS(dT + ddT))*ddT/2.0;
}

                // 2.9.4) Fatia de área obtida da função normalnucV pelo Método dos Trapézios

double trapezionormalnucV(double dT){   // Integral pelo método dos trapézios, de dn/d(dT) para nucleação no volume
        return  (normalnucV(dT) + normalnucV(dT + ddT))*ddT/2.0;
}


void NucleationStart(){


        #if dTvsigma != 0.0 || dTssigma != 0.0
        double dTi;
        double sumnuc; // Soma de núcleos
        double dTo;
        int startnuc=0 ; // "Switch que indica se o dT é inicial
        int nteste=0;
        int alea;

        #endif

        #ifdef NUCALLWALLDOWN
                nsmax = mx;
        #endif

        double dTn;

        int alea1, alea2;



        contadorS = 0;
        contadorV = 0;
        contadordraw = 0;

        contadorV = 0;
        contador = 0;

        #ifdef NUCLEATIONCENTER

        #if dTvsigma==0.0
        dTn = dTvnuc;
        printf("\nDistribuindo nucleos no interior do molde...");
        while( contadorV < mx*my && contadorV < (int) nvmax && contador < 10*mx*my){

                // Nucleação no centro

                alea1 = random(mx - NUCLEFT - NUCRIGHT) + NUCLEFT;
                alea2 = random(my - NUCUP - NUCDOWN) + NUCUP;
                //for(alea2=0;alea2<my;alea2++){
                //for(alea1=0;alea1<mx;alea1++){
                if(CA[(alea2)*mx + (alea1)].state == 0 && CA[(alea2)*mx + (alea1)].dT < 0.0){
                        CA[(alea2)*mx + (alea1)].dT = dTn;
                        contadorV++;
                }
                contador++;
                //}
                //
        }

        #ifdef NUCLEATIONCENTER2
        dTn = dTvnuc2;
        printf("\nDistribuindo 2o nucleos no interior do molde...");
        while( contadorV < mx*my && contadorV < (int) nvmax2 + (int) nvmax && contador < 10*mx*my){

                // Nucleação no centro

                alea1 = random(mx - NUCLEFT - NUCRIGHT) + NUCLEFT;
                alea2 = random(my - NUCUP - NUCDOWN) + NUCUP;
                //for(alea2=0;alea2<my;alea2++){
                //for(alea1=0;alea1<mx;alea1++){
                if(CA[(alea2)*mx + (alea1)].state == 0 && CA[(alea2)*mx + (alea1)].dT < 0.0){
                        CA[(alea2)*mx + (alea1)].dT = dTn;
                        contadorV++;
                }
                contador++;
                //}
                //
        }
        #endif
        #endif


        ia=ia;

        #if dTvsigma != 0.0

        sumnuc = 0;
        printf("\nDistribuindo nucleos no interior do molde...");
        for(dTi=0.0; dTi< dTvnuc + 300.0*dTvsigma; dTi+=(double) ddT){
                if(startnuc==0){
                        startnuc=1;
                        dTo=dTi;
                }
                sumnuc+=trapezionormalnucV(dTi);
                if(sumnuc>=1){


                        // Estratégia D3

                        if(dTo < dTvnuc - 5*dTvsigma){
                                dTo = dTvnuc - 5*dTvsigma;
                        }


                        nteste=(int) sumnuc;
                        dTn=(dTi + dTo)/2.0;    // Superresfriamento médio



                        srand(rand());

                        for(int k=1; k<=nteste && contadorV < (mx - NUCLEFT - NUCRIGHT)*(my - NUCUP - NUCDOWN); k++){              // Nucleação no centro
                                alea1 = random(mx - NUCLEFT - NUCRIGHT) + NUCLEFT;
                                alea2 = random(my - NUCUP - NUCDOWN) + NUCUP;
                                if(CA[(alea2)*mx + (alea1)].state == 0 && CA[(alea2)*mx + (alea1)].dT < 0.0){
                                        CA[(alea2)*mx + (alea1)].dT = dTn;
                                        contadorV++;
                                }
                                else{
                                        k=k-1;
                                }
                        }


                        startnuc=0;
                        sumnuc=0;
                }
        }
        #endif

        #endif



        
        #ifdef NUCLEATIONEqualNinVFs

        #if dTvsigma==0.0
        dTn = dTvnuc;
        printf("\nDistribuindo nucleos no interior do molde...");
        for(ia=0;ia<mmay;ia++){
                for(ja=0;ja<mmax;ja++){
                                contador = 0;
                                contadorV = 0;
                                while( contadorV < mmix*mmiy && contadorV < (int) nvmax && contador < 10*mmix*mmiy ){

                                                // Nucleação no centro
                                                int ipermitido = mmiy;
                                                int idesloc = 0;
                                                int jpermitido = mmix;
                                                int jdesloc = 0;

                                                if(ia == 0){
                                                                ipermitido -= NUCUP;
                                                                idesloc += NUCUP;
                                                }
                                                if(ia == mmay - 1){
                                                                ipermitido -= NUCDOWN;
                                                                idesloc -= NUCDOWN;
                                                }
                                                if(ja == 0){
                                                                jpermitido -= NUCLEFT;
                                                                jdesloc += NUCLEFT;
                                                }
                                                if(ja == mmax - 1){
                                                                jpermitido -= NUCRIGHT;
                                                                jdesloc -= NUCRIGHT;
                                                }
                                                alea1 = ja*mmix + random(jpermitido) + jdesloc;
                                                alea2 = ia*mmiy + random(ipermitido) + idesloc;

                                                if(CA[(alea2)*mx + (alea1)].state == 0 && CA[(alea2)*mx + (alea1)].dT < 0.0){
                                                CA[(alea2)*mx + (alea1)].dT = dTn;
                                                contadorV++;
                                                }
                                                contador++;

                                }
                }
        }
        #endif


        ia=ia;

        #if dTvsigma != 0.0

        sumnuc = 0;
        printf("\nDistribuindo nucleos no interior do molde...");
        for(dTi=0.0; dTi< dTvnuc + 300.0*dTvsigma; dTi+=(double) ddT){
                if(startnuc==0){
                        startnuc=1;
                        dTo=dTi;
                }
                sumnuc+=trapezionormalnucV(dTi);
                if(sumnuc>=1){


                        // Estratégia D3

                        if(dTo < dTvnuc - 5*dTvsigma){
                                dTo = dTvnuc - 5*dTvsigma;
                        }


                        nteste=(int) sumnuc;
                        dTn=(dTi + dTo)/2.0;    // Superresfriamento médio



                        srand(rand());

                        for(int k=1; k<=nteste; k++){              // Nucleação no centro
                                alea1 = random(mx - NUCLEFT - NUCRIGHT) + NUCLEFT;
                                alea2 = random(my - NUCUP - NUCDOWN) + NUCUP;
                                if(CA[(alea2)*mx + (alea1)].state == 0 && CA[(alea2)*mx + (alea1)].dT < 0.0){
                                        CA[(alea2)*mx + (alea1)].dT = dTn;
                                        contadorV++;
                                }
                                else{
                                        k=k-1;
                                }
                        }


                        startnuc=0;
                        sumnuc=0;
                }
        }
        #endif

        #endif




        #ifdef NUCLEATIONWALL
        #if dTssigma == 0.0
        #ifndef ALLNUCWALL
        dTn = dTsnuc;
        printf("\nDistribuindo nucleos na superficie do molde...");
        while( contadorS < mx && contadorS < nsmax && contador < 10*mx){

                // Nucleação no centro
                alea1 = random(mx);
                alea2 = random(my);
                if(CA[(my -1)*mx + (alea1)].state == 0 && CA[(my -1)*mx + (alea1)].dT < 0.0){
                        CA[(my -1)*mx + (alea1)].dT = dTn;
                        contadorS++;
                }
                contador++;
        }
        #endif
        #endif
        #endif

        #ifdef NUCLEATIONWALL
        #if dTssigma == 0.0
        #ifdef ALLNUCWALL
        dTn = dTvnuc;
        printf("\nDistribuindo nucleos na superficie do molde...");
        while( contadorV < mx*my && contadorV < nvmax && contador < 10*mx*my){

                // Nucleação no centro
                alea1 = random(mx);
                alea2 = random(my);
                if(volmicro[(alea2)*mx + (alea1)].state == 0 && volmicro[(alea2)*mx + (alea1)].dT < 0.0){
                        volmicro[(alea2)*mx + (alea1)].dT = dTn;
                        contadorV++;
                }
                contador++;
        }
        #endif
        #endif
        #endif

        #ifdef NUCLEATIONWALL
        #if dTssigma != 0.0
        sumnuc = 0;
        printf("\nDistribuindo nucleos na superficie do molde...");
        for(dTi=0.0; dTi< dTsnuc + 50.0*dTssigma ; dTi+=(double) ddT){ // Critério arbitrário
                if(startnuc==0){
                        startnuc=1;
                        dTo=dTi;
                }
                sumnuc+=trapezionormalnucS(dTi);

                if(sumnuc>=1){  // Se a integral naquela partição for maior que uma unidade...


                        // Estratégia D3

                        nteste=(int) sumnuc;
                        dTn=(dTi + dTo)/2.0;        // dT médio, após completar uma unidade na integral



                        srand(rand());
                        randomize();

                        #if NUCUP
                        int contadorup = 0;
                        for( int k=1; k<=nteste && contadorup <= 100*mx ; k++){     //Nucleação na parede superior


                                alea = random(mx - 1);

                                if(volmicro[0*mx + (alea + 1)].state == 0 && volmicro[0*mx + (alea + 1)].dT < 0.0){
                                        volmicro[0*mx + (alea + 1)].dT = dTn;
                                        contadorS = contadorS + 1;
                                        contadorup++ ;
                                }
                                else{
                                        k=k-1;
                                        contadorup++;
                                }
                        }
                        #endif

                        #if NUCRIGHT
                        int contadorright = 0;
                        for( int k=1; k<=nteste && contadorright <= 100*my ; k++){     //Nucleação na parede direita
                                alea = random(my - 1);
                                if(volmicro[(alea + 1)*mx + (mx - 1)].state == 0 && volmicro[(alea + 1)*mx + (mx - 1)].dT < 0.0){
                                        volmicro[(alea + 1)*mx + (mx - 1)].dT = dTn;
                                        contadorS = contadorS + 1;
                                        contadorright++;
                                }
                                else{
                                        k=k-1;
                                        contadorright++;
                                }
                        }
                        #endif

                        #if NUCDOWN
                        int contadordown = 0;
                        for( int k=1; k<=nteste && contadordown <= 100*mx; k++){     //Nucleação na parede inferior
                                alea = random(mx - 1);
                                if(CA[(my - 1)*mx + (mx - 1) - (alea + 1)].state == 0.0 && CA[(my - 1)*mx + (mx - 1) - (alea + 1)].dT < 0.0){
                                        CA[(my - 1)*mx + (mx - 1) - (alea + 1)].dT = dTn;
                                        contadorS = contadorS + 1;
                                        contadordown++;
                                }
                                else{
                                        k=k-1;
                                        contadordown++;
                                }
                        }
                        #endif

                        #if NUCLEFT
                        int contadorleft = 0;
                        for( int k=1; k<=nteste && contadorleft <= 100*my; k++){     // Nucleação na parede esquerda
                                alea = random(my - 1);
                                if(volmicro[((my - 1) - (alea + 1))*mx + 0].state == 0 && volmicro[((my - 1) - (alea + 1))*mx + 0].dT < 0.0){
                                        volmicro[((my - 1) - (alea + 1))*mx + 0].dT = dTn;
                                        contadorS = contadorS + 1;
                                        contadorleft++;
                                }
                                else{
                                        k=k-1;
                                        contadorleft++;
                                }
                        }
                        #endif

                        startnuc=0;
                        sumnuc=0;
                }
        }
        #endif
        #endif

        
        #ifdef NUCLEATIONEXACT

        for( int k = 0; k < SubstractsNumber ; k++){


                alea1 = random(mx - NUCLEFT - NUCRIGHT) + NUCLEFT;
                alea2 = random(my - NUCUP - NUCDOWN) + NUCUP;
                while( CA[alea2*mx + alea1].dT != -1.0 ){
                        alea1 = random(mx - NUCLEFT - NUCRIGHT) + NUCLEFT;
                        alea2 = random(my - NUCUP - NUCDOWN) + NUCUP;

                }

                CA[alea2*mx + alea1].dT = dTinstant;




        }
        #endif



}
