//---------------------------------------------------------------------------
#ifndef OutputFilesH
#define OutputFilesH
//---------------------------------------------------------------------------



#ifdef THERMALPROFILE
extern FILE *TMacro;
#endif

#ifdef SOLIDFRACTIONPROFILE
extern FILE *fSMacro;
#endif

#ifdef Eutetic
#ifdef EUTETICFRACTIONPROFILE
extern FILE *fEMacro;
#endif
#ifdef SUMSOLIDEUTETICFRACTIONPROFILE
extern FILE *fSfEMacro;
#endif
#endif

#ifdef THERMALHOPPROFILE
extern FILE *TMacroHop;
#endif

#ifdef SOLIDFRACTIONHOPPROFILE
extern FILE *fSMacroHop;
#endif

#ifdef PRINTGRAINS
extern FILE *gFile;
#endif

#ifdef PRINTCONTORNO
extern FILE *cFile;
#endif


#ifdef PRINTSOLUTEPROFILE
extern FILE *flClFile;
#endif

#ifdef PRINTLIQUIDFRACTIONPROFILE
extern FILE *fLMacro;
extern FILE * fLCAMacro;
#endif

#ifdef PRINTAverageVelocities
extern FILE *VLMacro;
extern FILE *VdMacro;
extern FILE *VsMacro;
extern FILE *elVLMacro;
#endif

#ifdef PRINTSURFACECD
extern FILE *SurfaceCd;
#endif

#ifdef PRINTVERTICALVALUES
extern FILE *VertCL;
extern FILE *VertVL;
extern FILE *VertVd;
extern FILE *VertVs;
extern FILE *VertfS;
extern FILE *VertflVl;
extern FILE *VertfsVs;
extern FILE *VertfdVd;
extern FILE *VertSe;
extern FILE *VertT;
extern FILE *VertCs;
extern FILE *VertC;
#endif

extern FILE *efvf;
extern FILE *elvl;

double AverageCsCdCl(double CVF[mmacrox*mmacroy][5]);

void InitiateColors(float red[nclassestheta], float green[nclassestheta], float blue[nclassestheta]);



#if PRINTGRAIN
extern FILE *gFile;
#endif



#ifdef THERMALPROFILE
void PrintTMacro();

#endif




#ifdef SOLIDFRACTIONPROFILE
void PrintfSMacro();

#endif

#ifdef Eutetic
#ifdef EUTETICFRACTIONPROFILE
void PrintfEMacro();
#endif
#ifdef SUMSOLIDEUTETICFRACTIONPROFILE
void PrintfSfEMacro();
#endif
#endif

#ifdef PRINTLIQUIDFRACTIONPROFILE
void PrintfLMacro();
#endif

#ifdef PRINTSOLUTEPROFILE


void PrintSoluteMacro();

void PrintSoluteMacro();
#endif

#ifdef THERMALHOPPROFILE
void PrintTMacroHop();
#endif

#ifdef SOLIDFRACTIONHOPPROFILE
void PrintfSMacroHop();
#endif

#ifdef PRINTAverageVelocities
void PrintAverageVelocities();
#endif



#ifdef PRINTCALL
void PrintCall(double CVF[mmacroy*mmacrox][5]);
#endif

#ifdef PRINTGRAINS

void PrintGrains();

#endif

#ifdef PRINTCONTORNO


void PrintContorno();

#endif



void Csolidaverage(double CsVF[mmacroy]);

void Caverage(double CVF[mmacroy]);


double CalculateSe( double La );


double Calculatefg(double Aa);


void CountLaAa(double *La, double *Aa);



void InitiatePrintCADE();

void AverageCdforVF();

void PrintCADE();





void PrintCADEComparisonResults();

void PrintProfiles();


void InitialVisual();


void InitiateFiles();

#ifdef PRINTVERTICALVALUES
void PrintVerticalValues(double CVF[mmacroy*mmacrox][5]);
#endif

#ifdef PRINTSURFACECD

void PrintSurfaceCd();

#endif
#endif
