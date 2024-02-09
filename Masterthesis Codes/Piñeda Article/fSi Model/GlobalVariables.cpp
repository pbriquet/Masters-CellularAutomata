//---------------------------------------------------------------------------
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#pragma hdrstop

//#include "GlobalVariables.h"
//---------------------------------------------------------------------------
#pragma package(smart_init)


#include "GlobalDefinitions.h"

#ifndef GlobalVariablesCPP
#define GlobalVariablesCPP


double T4vst[28][2];

#ifdef ThermoCouples
struct ThermoCoupleData{
        double Curvas_Resf[NUMBER_INICIAL_T][N_MEDIDAS];
        double TimeCurves[N_MEDIDAS];
        double PosInicial[NUMBER_INICIAL_T];//Posi��o dos termopares
};

struct ThermoCoupleData ThermoData;
#endif

double Vref;
bool ErroCAMotion;

double TCAdt;


struct ClusterMatrix{
        int thetaclass;     // A que classe de �ngulo theta o gr�o pertence (para manter a cor)
        int cluster;    // Indica a qual cluster o envelope pertence. De in�cio, � o mesmo que a classe. Ao formar um cluster, o �ndice do cluster � igual a classe menor dos envelopes envolvidos
        double V[2];    // Vetor velocidade Vx = V[0], Vy = V[1]
        double d[2];    // Vetor deslocamento dx = d[0], dy = d[1]
        bool left;      // Cluster deve se mover para a esquerda ? TRUE/FALSE
        bool right;     // Cluster deve se mover para a direita ? TRUE/FALSE
        bool up;        // Cluster deve se mover para cima ? TRUE/FALSE
        bool down;      // Cluster deve se mover para baixa ? TRUE/FALSE
        bool fixed;      // Cluster deve se tornar fixo ? TRUE/FALSE
        int size;       // C�lulas ocupadas pelo CLUSTER
        double fS;
        double Dave;

        #ifdef esiCluster
                double esi;
                double esidt;
                double dCddt;
                double Cdt;
                double de;
                double ve;
                int countint;
                double Ae;
                double Ve;
        #endif


};


struct MacroMatrix{

                double T;
                double dT;      // Pode n�o ser necess�rio ou evita a constru��o de uma matrix t + dt
                double dH;
                double CL;
                double CLdt;
                double VL;
                
                double Vd;
                double Vs;
                double Vf;
                double Se;

                double Cd;
                double fS;
                double dfS;
                #ifdef Eutetic
                double fE;
                double dfE;
                #endif
                double fL;
                double dfL;

                double fLCA;

                bool state;

                int Ne;
                int Nuc;



};


        // 1.4) Estrutura da matriz microsc�pica


struct MicroMatrix{

        double theta;   // Orienta��o da c�lula
        double L[4];       // Tamanho do semi-lado do envelope quadrado
        double dC[2];   // Vetor descentraliza��o do envelope
        int state;      // �ndice que indica o Estado f�sico da c�lula   0 = l�quido, 1 = molde, >=2 = �ndice do s�lido em crescimento, <-2 = �ndice do s�lido est�tico
        double fS; // Fra��o de s�lido
        double dfS;
        double Cd;
        double Cddt;
        double Cs;
        double Csdt;
        #ifdef eglocal
        double Amin;
        double Amax;
        double dfG;
        #endif

        #ifdef Eutetic
        double fE;
        double dfE;
        #endif

        #ifdef KeepClinCA
        double Cl;
        #endif

        double dT;

};




// Initialize Structures
struct MacroMatrix *VF;  /* Malha do modelo macrosc�pico */
void *auxpointmacro;

struct MacroMatrix *VFaux;      /* Malha do modelo macrosc�pico auxiliar */
void *auxpointmacroaux;

struct ClusterMatrix *clusters;      // O �ndica no registro corresponde ao "state" da malha - 1
void *auxpointclusters;

struct MicroMatrix *CA;  // Malha do modelo microsc�pico
void *auxpointmicro;

struct MicroMatrix *CAaux;  // Malha do modelo microsc�pico
void *auxpointmicroaux;

int NeKeeper[mmicrox*mmicroy*mmacrox*mmacroy];

double efvfyn[mmacrox*mmacroy];
double elvlyn[mmacrox*mmacroy];

double Se, eg;
struct MacroMatrix *VFP;
struct MacroMatrix *VFN;
struct MacroMatrix *VFS;
struct MicroMatrix *CAP;
struct MicroMatrix *CAN;
struct MicroMatrix *CAS;
struct MicroMatrix *CAE;
struct MicroMatrix *CAW;


int mmax = mmacrox;      // Quantidade de VF em uma linha
int mmay = mmacroy;      // Quantidade de VF em uma coluna

int mmix = mmicrox;      // Quantidade de CAs em uma linha dentro de um VF
int mmiy = mmicroy;      // Quantidade de CAs em uma coluna dentro de um VF

int mx = mmax*mmix;       // Quantidade de CAs em uma linha
int my = mmay*mmiy;       // Quantidade de CAs em uma coluna

double Lx = Lmalhax, Ly = Lmalhay;      // Tamanho do dom�nio em metros

double dxVF = Lx/mmax;        // Dimens�o horizontal de um VF
double dyVF = Ly/mmay;        // Dimens�o vertical de um VF

int NTotalGrains;

double dxCA = Lx/mx;    // Dimens�o horizontal de um CA
double dyCA = Ly/my;    // Dimens�o vertical de um VF

double Vo = 0.0;
double fSo = 0.0;

double t = 0.0;
double dt = dtmax;

double TVF;

double TCA;
double dHCA;
double CLVFDiffusion;
double CLVFTipGrowth;

double Rdl, Rp, CdlCA;
double fdP, velocity;

double fGCA;
double fGCAdt;


double gamaS1CA;
double aCA,bCA,cCA;

double gamaS2CA;
double gamaSCA;

double fSCA;
double CdCA;
#ifdef Eutetic
double fECA;
#endif

double JDdlCA;
double JDdlVF;
double CdJdlVF;
double JdlVF;
double JdlCA;

#ifdef CdJldExtraTerm
long double CdJldVl;
long double JldVl;
long double CdJldVs;
long double JldVs;

double AldOpen;
double AldOpenUp;
double AldOpenDown;
#endif


double CdAve;
double CdInterfaceAve;
double Cdmin;
double Cdmax;

double CsAve;

double dfSVF;
#ifdef Eutetic
double dfEVF;
#endif


double del_dl;

int phaseinterface = 0;


int stateholder = 0;
bool newgrain;

double vT[1+2*dimmax];
double aT[1+2*dimmax];
double vflCl[1+2*dimmax];
double aflCl[1 + 2*dimmax];
double vfl[1+2*dimmax];
double afl[1 + 2*dimmax];
double SUMaflvfl;
double SUMaflClvflCl;

double vfdCd[1+2*dimmax];
double afdCd[1 + 2*dimmax];
double vfd[1+2*dimmax];
double afd[1 + 2*dimmax];
double afdgamaS;
double aTgamaS;
double afdCdgamaS;
double SUMafdvfd;
double SUMafdCdvfdCd;
#ifdef Eutetic
double afdCdgamaE;
#endif

int side;

FILE *gamaSMacro;
FILE *gamaEMacro;

int ia, ja;
int i, j;
int ii, ji;


int contadorV;
int contadorS;
int contador;
int contadordraw;

int initialii = ((int) (my/2));
int initialji = (int) (mx/2);
int initialposition = ( initialii*mx) + initialji;


//double nvmax = pow( 6.0*pow(Nvmax,2.0)/M_PI, 1.0/3.0)*((mx - 2)*dxCA*(my - 2)*dyCA); //500;//267300*Lmalhax*Lmalhay; // //       /* �rea do "corte" analisado multiplicado pela densidade superficial equivalente (da volum�trica fornecida) de n�cleos */
//double nvmax = pow( 6.0*pow(Nvmax,2.0)/M_PI, 1.0/3.0)*Lx*Ly;
double nvmax = 6.5e4*Lx*Ly;
#ifdef NUCLEATIONCENTER2
double nvmax2 = 6.5e4*Lx*Ly;
#endif

//double nsmax = 583.0*Lx;     /* Per�metro do molde multiplicado pela densidade linear equivalente (da superficial fornecida) de n�cleos */
double nsmax = 2.0*sqrt(Nsmax/M_PI);

// float red[nclassestheta], green[nclassestheta], blue[nclassestheta]; // Color for each class of CA

double dtheta = M_PI/4.0/48.0;  // Angle discrete interval for CA classes

double A1, A2, Ar, Ah;






// Flow Partition Coefficient Calculation

#ifdef MartoranoLiquidBoundaryCondition

struct BoundaryConditionHolderStruct{
        double edVd;
        bool on;
        double vcluster;
};

struct BoundaryConditionHolderStruct BoundaryConditionHolder[mmicrox*mmacrox];
#endif

double DelAdjustmentValue;




double Kl(){


                #ifdef AlSi

                                #ifdef Co5
                                                return 176.0;
                                #endif

                                #ifdef Co7
                                                //return 36.5 + 0.028*TVF;
                                                return 60.5;
                                #endif
                #endif

                #ifdef AlCu

                                #ifdef Co3
                                                return 77.0;
                                #endif

                                #ifdef Co1
                                                return 90.0;
                                #endif

                                #ifdef Co10
                                                return 82.0;
                                #endif


                #endif


}

double Ks(){

                #ifdef AlSi

                                #ifdef Co5
                                                return 200.0;
                                #endif

                                #ifdef Co7
                                                //return 233.0 - 0.110*TVF;
                                                return 137.5;
                                #endif
                #endif

                #ifdef AlCu

                                #ifdef Co3
                                                return 153.0;
                                #endif

                                #ifdef Co1
                                                return 180.0;
                                #endif


                                #ifdef Co10
                                                return 160.0;
                                #endif


                #endif


}

#endif
