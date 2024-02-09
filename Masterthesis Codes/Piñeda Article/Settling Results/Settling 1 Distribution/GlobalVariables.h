//---------------------------------------------------------------------------
#ifndef GlobalVariablesH
#define GlobalVariablesH
//---------------------------------------------------------------------------

#include <stdio.h>

#include "GlobalDefinitions.h"

extern double T4vst[28][2];

extern double Vref;
extern bool ErroCAMotion;

#ifdef ThermoCouples
struct ThermoCoupleData{
        double Curvas_Resf[NUMBER_INICIAL_T][N_MEDIDAS];
        double TimeCurves[N_MEDIDAS];
        double PosInicial[NUMBER_INICIAL_T];//Posição dos termopares
};

extern struct ThermoCoupleData ThermoData;
#endif

struct ClusterMatrix{
        int thetaclass;     // A que classe de ângulo theta o grão pertence (para manter a cor)
        int cluster;    // Indica a qual cluster o envelope pertence. De início, é o mesmo que a classe. Ao formar um cluster, o índice do cluster é igual a classe menor dos envelopes envolvidos
        double V[2];    // Vetor velocidade Vx = V[0], Vy = V[1]
        double d[2];    // Vetor deslocamento dx = d[0], dy = d[1]
        bool left;      // Cluster deve se mover para a esquerda ? TRUE/FALSE
        bool right;     // Cluster deve se mover para a direita ? TRUE/FALSE
        bool up;        // Cluster deve se mover para cima ? TRUE/FALSE
        bool down;      // Cluster deve se mover para baixa ? TRUE/FALSE
        bool fixed;      // Cluster deve se tornar fixo ? TRUE/FALSE
        int size;       // Células ocupadas pelo CLUSTER
        double fS;
        double Dave;


};


struct MacroMatrix{

                double T;
                double dT;      // Pode não ser necessário ou evita a construção de uma matrix t + dt
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


        // 1.4) Estrutura da matriz microscópica


struct MicroMatrix{

        double theta;   // Orientação da célula
        double L[4];       // Tamanho do semi-lado do envelope quadrado
        double dC[2];   // Vetor descentralização do envelope
        int state;      // Índice que indica o Estado físico da célula   0 = líquido, 1 = molde, >=2 = índice do sólido em crescimento, <-2 = índice do sólido estático
        double fS; // Fração de sólido
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




double Kl();
double Ks();


// Initialize Structures
extern struct MacroMatrix *VF;  /* Malha do modelo macroscópico */
extern void *auxpointmacro;

extern struct MacroMatrix *VFaux;      /* Malha do modelo macroscópico auxiliar */
extern void *auxpointmacroaux;

extern struct ClusterMatrix *clusters;      // O índica no registro corresponde ao "state" da malha - 1
extern void *auxpointclusters;

extern struct MicroMatrix *CA;  // Malha do modelo microscópico
extern void *auxpointmicro;

extern struct MicroMatrix *CAaux;  // Malha do modelo microscópico
extern void *auxpointmicroaux;

extern int NeKeeper[mmicrox*mmicroy*mmacrox*mmacroy];

extern double efvfyn[mmacrox*mmacroy];
extern double elvlyn[mmacrox*mmacroy];

extern double Se, eg;
extern struct MacroMatrix *VFP;
extern struct MacroMatrix *VFN;
extern struct MacroMatrix *VFS;
extern struct MicroMatrix *CAP;
extern struct MicroMatrix *CAN;
extern struct MicroMatrix *CAS;
extern struct MicroMatrix *CAE;
extern struct MicroMatrix *CAW;


extern int mmax;      // Quantidade de VF em uma linha
extern int mmay;      // Quantidade de VF em uma coluna

extern int mmix;      // Quantidade de CAs em uma linha dentro de um VF
extern int mmiy;      // Quantidade de CAs em uma coluna dentro de um VF

extern int mx;       // Quantidade de CAs em uma linha
extern int my;       // Quantidade de CAs em uma coluna

extern double Lx, Ly;      // Tamanho do domínio em metros

extern double dxVF ;        // Dimensão horizontal de um VF
extern double dyVF ;        // Dimensão vertical de um VF

extern int NTotalGrains;

extern double dxCA ;    // Dimensão horizontal de um CA
extern double dyCA ;    // Dimensão vertical de um VF

extern double Vo;
extern double fSo;

extern double t;
extern double dt;

extern double TVF;

extern double TCA;
extern double dHCA;
extern double CLVFDiffusion;
extern double CLVFTipGrowth;

extern double Rdl, Rp, CdlCA;
extern double fdP, velocity;

extern double fGCA;
extern double fGCAdt;


extern double gamaS1CA;
extern double aCA,bCA,cCA;

extern double gamaS2CA;
extern double gamaSCA;

extern double fSCA;
extern double CdCA;
#ifdef Eutetic
extern double fECA;
#endif

extern double JDdlCA;
extern double JDdlVF;
extern double CdJdlVF;
extern double JdlVF;
extern double JdlCA;

#ifdef CdJldExtraTerm
extern long double CdJldVl;
extern long double JldVl;
extern long double CdJldVs;
extern long double JldVs;

extern double AldOpen;
extern double AldOpenUp;
extern double AldOpenDown;
#endif


extern double CdAve;
extern double CdInterfaceAve;
extern double Cdmin;
extern double Cdmax;

extern double CsAve;

extern double dfSVF;
#ifdef Eutetic
extern double dfEVF;
#endif


extern double del_dl;

extern int phaseinterface;


extern int stateholder;
extern bool newgrain;

extern double vT[1+2*dimmax];
extern double aT[1+2*dimmax];
extern double vflCl[1+2*dimmax];
extern double aflCl[1 + 2*dimmax];
extern double vfl[1+2*dimmax];
extern double afl[1 + 2*dimmax];
extern double SUMaflvfl;
extern double SUMaflClvflCl;

extern double vfdCd[1+2*dimmax];
extern double afdCd[1 + 2*dimmax];
extern double vfd[1+2*dimmax];
extern double afd[1 + 2*dimmax];
extern double afdgamaS;
extern double aTgamaS;
extern double afdCdgamaS;
extern double SUMafdvfd;
extern double SUMafdCdvfdCd;
#ifdef Eutetic
extern double afdCdgamaE;
#endif

extern int side;

extern FILE *gamaSMacro;
extern FILE *gamaEMacro;

extern int ia, ja;
extern int i, j;
extern int ii, ji;


extern int contadorV;
extern int contadorS;
extern int contador;
extern int contadordraw;

extern int initialii;
extern int initialji;
extern int initialposition;


//double nvmax = pow( 6.0*pow(Nvmax,2.0)/M_PI, 1.0/3.0)*((mx - 2)*dxCA*(my - 2)*dyCA); //500;//267300*Lmalhax*Lmalhay; // //       /* Área do "corte" analisado multiplicado pela densidade superficial equivalente (da volumétrica fornecida) de núcleos */
extern double nvmax;

extern double nsmax;     /* Perímetro do molde multiplicado pela densidade linear equivalente (da superficial fornecida) de núcleos */

// float red[nclassestheta], green[nclassestheta], blue[nclassestheta]; // Color for each class of CA

extern double dtheta;  // Angle discrete interval for CA classes

extern double A1, A2, Ar, Ah;






// Flow Partition Coefficient Calculation

#ifdef MartoranoLiquidBoundaryCondition

struct BoundaryConditionHolderStruct{
        double edVd;
        bool on;
        double vcluster;
};

extern struct BoundaryConditionHolderStruct BoundaryConditionHolder[mmicrox*mmacrox];
#endif

extern double DelAdjustmentValue;


#endif
