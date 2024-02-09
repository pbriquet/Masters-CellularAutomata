//---------------------------------------------------------------------------
#ifndef PhaseTransfH
#define PhaseTransfH
//---------------------------------------------------------------------------






void dHDiffusive();

void getCoefficientdfSMicro();

void getCoefficientsdCLMacro();


void dfSMacro();

void dTMacro();

void dCLMacro();

void dfLMacro();


#ifdef Eutetic
void dfSMicro();
#endif


#ifndef Eutetic
void dfSMicro();

#endif

#endif
