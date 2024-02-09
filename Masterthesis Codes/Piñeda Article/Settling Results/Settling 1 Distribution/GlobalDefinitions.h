//---------------------------------------------------------------------------
#ifndef GlobalDefinitionsH
#define GlobalDefinitionsH
//---------------------------------------------------------------------------



// Switchs

        // PreSetup


                #define NoMaterialPreSetup
                #define NoMeshPreSetup
                #define NoInitialPreSetup
                #define NoNucleationPreSetup
                #define NoSubstractPreSetup
                #define NoHeatTransferPreSetup
                #define NoAlgorithmPreSetup
                #define NoLiquidDiffusionPreSetup
                #define NoTipGrowthPreSetup
                #define NoCACapturePropertiesPreSetup
                #define NoMotionPreSetup
                #define NoResultsPreSetup

                                //#define DeterministicModel
                                //#define CADEModel
                                //#define GandinModel
                                //#define ModeltoGandin

                #define CADEArticlePreSetup
                //#define ViniciusTFPreSetup
                //#define WangBeckermann1994PreSetup

                #define GrowthAlgorithm
                #define NucleationAlgorithm

                #include "SetupDatabase.h"
                

        #ifdef NoMaterialPreSetup
        // Material
                #define AlSi
        // Composition
                #define Co7
        #endif

        // Mesh and Domain Setup
                #ifdef NoMeshPreSetup

                                #define mmacrox 1
                                #define mmacroy 32
                                #define mmicrox 66
                                #define mmicroy 9

                                #define Lmalhax 0.0275
                                #define Lmalhay 0.12

                                #define dimmax 3
                                #define dimCA 2
                                #define dimVF 1
                #endif

        // Algorithm

                #ifdef NoTipGrowthPreSetup
                        //#define SOLUTALVELOCITY
                        #define THERMALVELOCITY
                        #define DendriticTipEquilibriumWithVFTemperature
                        #define Cofarfield
                        //#define CLVFTipGrowthfarfield
                        //#define CLInterpolfarfield
                        //#define CONSTANTVELOCITY
                        //#define THERMALVELOCITY2
                        #define THERMALVELOCITY3
                #endif

                #ifdef NoLiquidDiffusionPreSetup
                        #define CLinCAsVFforLiquidDiffusion
                        //#define CLInterpolforLiquidDiffusion
                        //#define CLequalCoforLiquidDiffusion

                        //#define CLVFDiffusionDiffusion
                        //#define BlockinterVFDiffusion
                        //#define DiffusionThroughVFs
                        //#define DiffusionDistanceColumnar
                        //#define inWallThereisInterExtraDiffusion

                        //#define STEREOMICRO
                        //#define STEREOMACRO
                        //#define STEREOMACROGROWTH
                        //#define STEREOOneSphereMACRO
                        //#define STEREOOneOctahedronMACRO


                        //#define SeCACorrection
                        //#define SeCACorrectionDiagonal
                        //#define SeCACorrectionDiagonalSolutalDiffusionLenght
                        //#define SeCACorrectionSphere
                        #define NSidesInLiquidDiffusion
                        //#define NSidesInGrainGrowth

                        #define InterliquidDiffusion DiffLiq
                        #define ExtraliquidDiffusion DiffLiq
                        //#define InterliquidDiffusion 0.0
                        //#define ExtraliquidDiffusion 0.0
                #endif

                #ifdef NoCACapturePropertiesPreSetup
                        //#define Cdl CA[ii*mx + ji].Cd
                        //#define Cld CA[ii*mx + ji].Cd
                        //#define Cdl (VF[ia*mmax + ja].T - Tf)/ml
                        //#define Cld (VF[ia*mmax + ja].T - Tf)/ml
                        //#define Cdl (TCA - Tf)/ml
                        //#define Cld (TCA - Tf)/ml
                        //#define Cdl VF[ia*mmax + ja].CL
                        //#define Cld VF[ia*mmax + ja].CL
                        #define Cld CdlCA
                        #define Cdl CdlCA

                        //#ifdef Cld==CdlCA
                        #define CdlEqualFluxinDiffusion
                        #define TEqualFluxinDiffusion
                        //#endif

                #endif

                #ifdef NoAlgorithmPreSetup

                        #define Eutetic

                        //#define NOGrowth
                        #define NOSolid

                        #define ConservativeLocalSolidfraction  // Model Solidfraction Calculus
                        //#define NonConservativeLocalSolidfraction
                        //#define GandinSolidfraction

                        //#define CddtThermalCAEquilibriumOnlyInFuture
                        //#define CddtThermalCAEquilibriumAlways
                        #define CddtConservative
                        //#define CddtNonConservative
                        //#define CddtEquilibriumVFt
                        //#define StorageTinCA
                        //#define TinEquilibriuminPhaseTransformation
                        #define TInterpolinPhaseTransformation

                        //#define CLdtConservative
                        #define CLdtSemiConservative
                        //#define CLdtNonConservative

                        #define CASolidConcentration

                        //#define NoInterpol

                        #define MotionAlgorithm
                                //#define NOMove

                                #ifdef MotionAlgorithm

                                #define MotionStabilityVerification
                                #define LevenspielVelocity
                                //#define ConstantTerminalVelocity
                                //#define AhujaVelocity

                                //#define SphereCA
                                #define OctahedronCA

                                #define MacrossegragetationCalculation
                                        #ifdef MacrossegragetationCalculation
                                                //#define FlowPartitionCoefficient
                                                #ifdef FlowPartitionCoefficient
                                                                #define efVfPequalseSvSP // (efvf)P=-(esvs)P OR
                                                                //#define efVfPequalseSvSN // (efvf)P=-(esvs)N OR

                                                                //#define lambda2 10.0*1.0e-6
                                                #endif
                                                #define MartoranoLiquidBoundaryCondition
                                                        #ifdef MartoranoLiquidBoundaryCondition
                                                                #define HoldBoundaryCondition
                                                        #endif
                                                #define CdJldExtraTerm
                                                #define CLConvectiveTermCalculation
                                                                #ifdef CLConvectiveTermCalculation
                                                                                #define fLCLupwind
                                                                                #define vLupwind
                                                                                //#define vLCentral

                                                                                #define RefreshfLdtvalueConvective
                                                                #endif
                                                //#define JldVsDiscrete
                                                //#define KeepClinCA

                                                // Which del_dl to use?

                                                //#define GrowthDel
                                                //#define ConvecDel
                                                //#define RanzMarshallDel
                                                #define NewDeldl

                                                //#define CdlCAConstante
                                                        #define CdlCAConst Co
                                                //#define CdCAConstante
                                                        #define CdCAConst Co
                                                //#define fSCAConstante
                                                        #define fSCAConst 0.0

                                                //#define NOPhaseTransformation
                                                //#define SolidSphere

                                                //#define sizeofInitialGrain dxCA
                                                #define sizeofInitialGrain Lx/4.0
                                                #define egproj 0.8
                                                //#define sizeofInitialGrain Lx*egproj/NumberofInitialGrainsx

                                                #define DelAdjustment
                                                        #ifdef DelAdjustment
                                                                #define DelAdjustConstant 1000.0
                                                        #define SelectiveDelAdjustment
                                                        #endif
                                        #endif

                                #endif


                        //#define eglocal

                #endif

                #ifdef NoMotionPreSetup
                        #define OCTAHEDRON
                        //#define SPHERE
                #endif


                #ifdef NoHeatTransferPreSetup

                                //#define HomogenousVFTemperature
                                //#define VFHeatBalance
                                //#define ZeroDimensions
                                //#define VirtualZeroDimensions
                                #define ThermoCouples

                                #define ScheilMicrofS

                                #define CoolingRate -20.0 // K/s

                                #define Tfarfieldup 298.0         // Distant Temperature
                                #define Tfarfielddown 298.0 //Tfar()         // Distant Temperature
                                #define Tfarfieldleft 298.0
                                #define Tfarfieldright 298.0
                                #define Hconvectiondown 250.0 //Hvar() //250.0       // H of convection
                                #define Hconvectionup 0.0 //250.0
                                #define Hconvectionleft 0.0
                                #define Hconvectionright 0.0
                #endif


                // Initial Conditions and Rum Time Setup
                #ifdef NoInitialPreSetup

                        //#define dTinitial ml*(CdCAConst - Co)
                        #define dTinitial 100.0
                        #define To Tliq + dTinitial

                        #define tinitial 0.0
                        #define dtmax 0.002
                        #define tmax 364.0

                        //#define tmaxperV

                #endif



                // Nucleation Setup
                #ifdef NoNucleationPreSetup

                         //#define LogNormalDistribution
                                #define ShapeParameter 0.7
                        #define NormalDistribution

                        #define Nsmax 3.0e6             // Number of substracts in mould boundary
                                #define dTsnuc  0.0             // Boundary Substracts Mean supercooling
                                #define dTssigma 0.0            // Boundary Substracts standard deviation of supercooling
                                #define dTvnuc 6.5              // Volume Substracts Mean supercooling
                                        #ifdef NormalDistribution
                                                #define dTvsigma 0.0          // Volume Substracts standard deviation of supercooling
                                        #endif
                                        #ifdef LogNormalDistribution
                                                //#define dTvsigma 0.05
                                        #endif
                                #define Nvmax 5.0e6            // Number of substracts in mould volume
                                #define dTinstant 6.5

                                #define SubstractsNumber 30
                #endif
                
                        #define ddT 0.00001             // Passo de super-resfriamento para integração numérica
                        #define dTsolidtotal 3456.0     // Maximum dT for solidification of the whole mould
                        #define nclassestheta 48
                        #define nucleosmax mx*my*(1 - 0.1)



                #ifdef NoSubstractPreSetup


                        //#define NUCLEATIONEXACT
                        #define NUCLEATIONCENTER
                                //#define NUCLEATIONEqualNinVFs
                                        #define NequalinVF 1
                        #define NUCLEATIONWALL
                        #define NUCLEFT 0
                        #define NUCRIGHT 0
                        #define NUCUP 0
                        #define NUCDOWN 0
                        //#define NUCALLWALLUP
                        #define NUCALLWALLDOWN

                        //#define INITIALGRAIN
                        //#define INITIALGRAINS
                                #define NumberofInitialGrainsx 5
                                #define NumberofInitialGrainsy 5

                        //#define ALL1GRAIN

                #endif

                // Results

                #ifdef NoResultsPreSetup
                //#define PRINTCADERESULTS
                #define CADEResultsTimeInterval dt
                #define PRINTGRAINS
                        // #define PRINTGrainsClusters
                //#define PRINTCONTORNO
                #define GrainsTimeInterval      10.0
                #define PRINTPROFILES
                #define ProfilesTimeInterval    0.1
                #define PRINTTIME
                #define TimeTimeInterval        dt

                                #define THERMALPROFILE
                                #define SOLIDFRACTIONPROFILE
                                #define EUTETICFRACTIONPROFILE
                                #define SUMSOLIDEUTETICFRACTIONPROFILE
                                #define PRINTSOLUTEPROFILE
                                //#define THERMALHOPPROFILE
                                //#define SOLIDFRACTIONHOPPROFILE
                                #define PRINTLIQUIDFRACTIONPROFILE
                                #define PRINTAverageVelocities

                                //#define PRINTSURFACECD

                                #define PRINTSeProfile
                                #define PRINTVERTICALVALUES
                                #define PRINTCALL
                #endif
                //#define PRINTHOPPROFILES

                #define ifProfilesTimeInterval (t >= 0.0 && t < 0.0 + dt) || t*(1.0/ProfilesTimeInterval) - (int) (t*(1.0/ProfilesTimeInterval) ) >= 0.0 && t*(1.0/ProfilesTimeInterval) - (int) (t*(1.0/ProfilesTimeInterval) ) < dt*(1.0/ProfilesTimeInterval)
                #define ifGrainsTimeInterval (t >= 0.0 && t < 0.0 + dt) || t*(1.0/GrainsTimeInterval) - (int) (t*(1.0/GrainsTimeInterval) ) >= 0.0 && t*(1.0/GrainsTimeInterval) - (int) (t*(1.0/GrainsTimeInterval) ) < dt*(1.0/GrainsTimeInterval)
                //#define ifGrainsTimeInterval (t >= 0.0 && t < 0.0 + dt) || (t >= 0.75 && t < 0.75 + dt) || (t >= 1.35 && t < 1.35 + dt) || (t >= 2.10 && t < 2.10 + dt)
                #define ifCADEResultsTimeInterval (t >= 0.0 && t < 0.0 + dt) || t*(1.0/CADEResultsTimeInterval) - (int) (t*(1.0/CADEResultsTimeInterval) ) >= 0.0 && t*(1.0/CADEResultsTimeInterval) - (int) (t*(1.0/CADEResultsTimeInterval) ) < dt*(1.0/CADEResultsTimeInterval)
                #define ifTimeTimeInterval (t >= 0.0 && t < 0.0 + dt) || t*(1.0/TimeTimeInterval) - (int) (t*(1.0/TimeTimeInterval) ) >= 0.0 && t*(1.0/TimeTimeInterval) - (int) (t*(1.0/TimeTimeInterval) ) < dt*(1.0/TimeTimeInterval)





#ifdef SPHERE
#define K1 0.7554
#define K2 0.8243
#endif

#ifdef OCTAHEDRON
#define K1 1.1272
#define K2 0.9677
#endif

// Properties Abbreviations
#define p soliddensity
#define Kliquid Kl()   // Conduction Parameter K for liquid phase (W/(m.K))
#define Ksolid Ks()    // Conduction Parameter K for solid phase (W/(m.K))
#define Dl DiffLiq



// General Properties
#define sigma 0.0253    // Stability Coefficient
#define gravity 9.8             // Gravity

#define fSmax 0.99999999
#define fdminVF 1.0/mmix/mmiy

#define fGmax 0.99999
#define Lmin 1e-10
//#define fLVFmin 2.0/mmicrox/mmicroy

#ifdef eglocal
#define fLVFmin 1.0e-8
#endif
#ifndef eglocal
#define fLVFmin 0.01/mmix/mmiy
#endif
                                             
#define TruncarfLdt
//#define TruncarfLVelocidade



// ................. definições para calculo de troca de calor ...........
#define Pos0 0.0     //(m)  posição 0
#define Pos1 0.024    //(m)  posição 1
#define Pos2 0.041    //(m)  posição 2
#define Pos3 0.059    //(m)  posição 3
#define Pos4 0.08    //(m)  posição 4
#define Pos5 0.099    //(m)  posição 5
#define Pos6 0.12   //0.121   //(m)  posição 6

#define NUMBER_INICIAL_T 7 //Number of cooling curves
#define N_MEDIDAS 3650  //Número de posições na curva de resfriamento













// Abbreviations
#define VfP VF[ia*mmax + ja]
#define VfN VF[(ia - 1)*mmax + ja]
#define VfS VF[(ia + 1)*mmax + ja]
#define VfE VF[ia*mmax + ja + 1]
#define VfW VF[ia*mmax + ja - 1]
#define CaP CA[ii*mx + ji]
#define CaN CA[(ii - 1)*mx + ji]
#define CaS CA[(ii + 1)*mx + ji]
#define CaE CA[ii*mx + ji + 1]
#define CaW CA[ii*mx + ji - 1]
#define RunVFLine               ia = 0; ia<mmay ; ia++
#define RunVFCol                ja = 0; ja<mmax ; ja++
#define RunCAinVFLine           ii = ia*mmiy; ii < (ia + 1)*mmiy; ii++
#define RunCAinVFCol            ji = ja*mmix; ji < (ja + 1)*mmix; ji++
#define RunCALine               ii = 0; ii < my ; ii++
#define RunCACol                ji = 0; ji < mx ; ji++
#define RunInverseVFLine               ia = mmay - 1; ia>=0 ; ia--
#define RunInverseVFCol                ja = mmax - 1; ja>=0 ; ja--
#define RunInverseCAinVFLine           ii = (ia + 1)*mmiy - 1; ii >= ia*mmiy ; ii--
#define RunInverseCAinVFCol            ji = (ja + 1)*mmix - 1; ji >= ja*mmiy ; ji--
#define RunInverseCALine               ii = my - 1; ii >= 0 ; ii--
#define RunInverseCACol                ji = mx - 1; ji >= 0 ; ji--
#define RunTime                 t = tinitial ; t <= tmax ; t+=dt
#define ClusterP clusters[ ((int) fabs( CaP.state )) - 2].cluster
#define ClusterE clusters[ ((int) fabs( CaE.state )) - 2].cluster
#define ClusterS clusters[ ((int) fabs( CaS.state )) - 2].cluster
#define ClusterN clusters[ ((int) fabs( CaN.state )) - 2].cluster
#define ClusterW clusters[ ((int) fabs( CaW.state )) - 2].cluster
#define RunClusterMatrix        ii = 0; ii < NTotalGrains; ii++

#define CaPaux CAaux[ii*mx + ji]
#define CaNaux CAaux[(ii - 1)*mx + ji]
#define CaSaux CAaux[(ii + 1)*mx + ji]
#define CaEaux CAaux[ii*mx + ji + 1]
#define CaWaux CAaux[ii*mx + ji - 1]

#define VfPaux VFaux[ia*mmax + ja]
#define VfNaux VFaux[(ia - 1)*mmax + ja]
#define VfSaux VFaux[(ia + 1)*mmax + ja]
#define VfEaux VFaux[ia*mmax + ja + 1]
#define VfWaux VFaux[ia*mmax + ja - 1]

#define ClusterPaux clusters[ ((int) fabs( CaPaux.state )) - 2].cluster

#define FuturePosition (ii - clusters[ClusterPaux].up + clusters[ClusterPaux].down)*mx + ji - clusters[ClusterPaux].left + clusters[ClusterPaux].right

#define FixClusterConditions ii == my - 1 || ii == 0 // || ji == 0 || ji == mx - 1
//#define FixClusterConditions 1 == 0


#define EAST 0
#define NORTH 1
#define WEST 2
#define SOUTH 3





#ifdef CADEArticlePreSetup

        #ifndef NoAlgorithmPreSetup

        #endif

        #ifndef NoMeshPreSetup
                #define mmacrox 1
                #define mmacroy 1
                #define mmicrox 200
                #define mmicroy 200

                #define Lmalhax 0.0079
                #define Lmalhay 0.0079

                #define dimmax 3
                #define dimCA 2
                #define dimVF 1
        #endif

        #ifndef NoHeatTransferPreSetup

                //#define HomogenousVFTemperature
                //#define VFHeatBalance
                #define ZeroDimensions
                //#define VirtualZeroDimensions

                #define CoolingRate -45.0 // K/s

                #define Tfarfieldup 298.0         // Distant Temperature
                #define Tfarfielddown 298.0         // Distant Temperature
                #define Tfarfieldleft 298.0
                #define Tfarfieldright 298.0
                #define Hconvectiondown 250.0 //250.0       // H of convection
                #define Hconvectionup 0.0 //250.0
                #define Hconvectionleft 0.0
                #define Hconvectionright 0.0

        #endif


        #ifndef NoInitialPreSetup
                #define dTinitial 0.0
                #define To Tliq + dTinitial
                #define tinitial 0.0
                #define dtmax 0.001
                #define tmax 8.0
        #endif

        #ifndef NoNucleationPreSetup
                #define Nsmax 2.5e8             // Number of substracts in mould boundary
                #define dTsnuc  0.5             // Boundary Substracts Mean supercooling
                #define dTssigma 0.1            // Boundary Substracts standard deviation of supercooling
                #define dTvnuc 0.0              // Volume Substracts Mean supercooling
                #define dTvsigma 0.0            // Volume Substracts standard deviation of supercooling
                #define Nvmax 5.5e8            // Number of substracts in mould volume

                #define dTinstant 0.0

                #define SubstractsNumber 30
        #endif

        #ifndef NoMaterialPreSetup
                #define AlSi
                #define Co 5.0
        #endif

        #ifndef NoSubstractPreSetup
                        #define NUCLEATIONEXACT
                        //#define NUCLEATIONCENTER
                        //#define NUCLEATIONWALL
                        #define NUCLEFT 0
                        #define NUCRIGHT 0
                        #define NUCUP 0
                        #define NUCDOWN 0
                        //#define NUCALLWALLDOWN

                        //#define INITIALGRAIN


        #endif

                // Results

        #ifndef NoResultsPreSetup
                #define PRINTCADERESULTS
                #define CADEResultsTimeInterval 0.01
                #define PRINTGRAINS
                #define PRINTCONTORNO
                #define GrainsTimeInterval      1.0
                #define PRINTPROFILES
                #define ProfilesTimeInterval    0.01
                #define PRINTTIME
                #define TimeTimeInterval        CADEResultsTimeInterval

                                #define THERMALPROFILE
                                #define SOLIDFRACTIONPROFILE
                                #define EUTETICFRACTIONPROFILE
                                #define PRINTSOLUTEPROFILE
                                //#define THERMALHOPPROFILE
                                //#define SOLIDFRACTIONHOPPROFILE
                                #define PRINTLIQUIDFRACTIONPROFILE



        #endif
#endif

#ifdef ViniciusTFPreSetup

                #ifndef NoAlgorithmPreSetup

                #endif

                #ifndef NoMeshPreSetup
                                #define mmacrox 1
                                #define mmacroy 20
                                #define mmicrox 110
                                #define mmicroy 11

                                #define Lmalhax 0.07
                                #define Lmalhay 0.14

                                #define dimmax 3
                                #define dimCA 2
                                #define dimVF 1
                #endif

                #ifndef NoNucleationPreSetup
                                #define Nsmax 2.5e8             // Number of substracts in mould boundary
                                #define dTsnuc  0.5             // Boundary Substracts Mean supercooling
                                #define dTssigma 0.1            // Boundary Substracts standard deviation of supercooling
                                #define dTvnuc 4.5              // Volume Substracts Mean supercooling
                                #define dTvsigma 0.5            // Volume Substracts standard deviation of supercooling
                                #define Nvmax pow( pow( 267300.0 , 3.0)*M_PI/6.0 , 0.5)            // Number of substracts in mould volume
                                #define dTinstant 0.0

                                #define SubstractsNumber 30
                #endif

                #ifndef NoMaterialPreSetup
                                #define AlSi
                                #define Co 7.0
                #endif

                #ifndef NoInitialPreSetup

                        #define dTinitial 100.0
                        #define To Tliq + dTinitial

                        #define tinitial 0.0
                        #define dtmax 0.002
                        #define tmax 1400.0

                #endif

                #ifndef NoSubstractPreSetup
                        //#define NUCLEATIONEXACT
                        #define NUCLEATIONCENTER
                        #define NUCLEATIONWALL
                        #define NUCLEFT 0
                        #define NUCRIGHT 0
                        #define NUCUP 0
                        #define NUCDOWN 1

                #endif

                #ifndef NoHeatTransferPreSetup

                                //#define HomogenousVFTemperature
                                #define VFHeatBalance
                                //#define ZeroDimensions
                                //#define VirtualZeroDimensions

                                #define CoolingRate -45.0 // K/s

                                #define Tfarfieldup 298.0         // Distant Temperature
                                #define Tfarfielddown 298.0         // Distant Temperature
                                #define Tfarfieldleft 298.0
                                #define Tfarfieldright 298.0
                                #define Hconvectiondown 250.0 //250.0       // H of convection
                                #define Hconvectionup 0.0 //250.0
                                #define Hconvectionleft 0.0
                                #define Hconvectionright 0.0

                #endif

                #ifndef NoResultsPreSetup
                                //#define PRINTCADERESULTS
                                #define CADEResultsTimeInterval 0.01
                                #define PRINTGRAINS
                                #define GrainsTimeInterval      50.0
                                #define PRINTPROFILES
                                #define ProfilesTimeInterval    1.0
                                #define PRINTTIME
                                #define TimeTimeInterval        0.1



                                #define THERMALPROFILE
                                #define SOLIDFRACTIONPROFILE
                                #define EUTETICFRACTIONPROFILE
                                #define PRINTSOLUTEPROFILE
                                #define THERMALHOPPROFILE
                                #define SOLIDFRACTIONHOPPROFILE

                #endif


#endif


#ifdef WangBeckermann1994PreSetup

                #ifndef NoMaterialPreSetup
                                #define AlCu
                                #define Co 3.0
                #endif

                #ifndef NoMeshPreSetup
                                #define mmacrox 1
                                #define mmacroy 20
                                #define mmicrox 110
                                #define mmicroy 11

                                #define Lmalhax 0.05
                                #define Lmalhay 0.10

                                #define dimmax 3
                                #define dimCA 2
                                #define dimVF 1
                #endif

                #ifndef NoHeatTransferPreSetup

                                //#define HomogenousVFTemperature
                                #define VFHeatBalance
                                //#define ZeroDimensions
                                //#define VirtualZeroDimensions

                                //#define CoolingRate -45.0 // K/s

                                #define Tfarfieldup 298.0         // Distant Temperature
                                #define Tfarfielddown Tfar()         // Distant Temperature
                                #define Tfarfieldleft 298.0
                                #define Tfarfieldright 298.0
                                #define Hconvectiondown 2000.0450 //250.0       // H of convection
                                #define Hconvectionup 0.0 //250.0
                                #define Hconvectionleft 0.0
                                #define Hconvectionright 0.0

                #endif

                #ifndef NoNucleationPreSetup
                                #define Nsmax 2.5e8             // Number of substracts in mould boundary
                                #define dTsnuc  0.0             // Boundary Substracts Mean supercooling
                                #define dTssigma 0.0            // Boundary Substracts standard deviation of supercooling
                                #define dTvnuc 0.0              // Volume Substracts Mean supercooling
                                #define dTvsigma 0.0            // Volume Substracts standard deviation of supercooling
                                #define Nvmax 1.0e5            // Number of substracts in mould volume
                                #define dTinstant 0.0

                                #define SubstractsNumber 30
                #endif

                #ifndef NoInitialPreSetup

                        #define dTinitial 50.0
                        #define To Tliq + dTinitial

                        #define tinitial 0.0
                        #define dtmax 0.002
                        #define tmax 4000.0

                #endif

                #ifndef NoSubstractPreSetup
                        //#define NUCLEATIONEXACT
                        #define NUCLEATIONCENTER
                        #define NUCLEATIONWALL
                        #define NUCLEFT 0
                        #define NUCRIGHT 0
                        #define NUCUP 0
                        #define NUCDOWN 1
                        #define NUCALLWALLDOWN

                #endif

#endif



#ifdef AlSi


        #ifdef Co5

                #define Co 5.0

                #define Tliq 894.5         // Liquidus Temperature for the given composition (7% w Si) (Kelvin)
                #define Tf 933.0
                #define Teut 850.0         // Eutetic Temperature (Estimated Value) (ºC)

                #define a2 0.00000290           // Estimated Parameter for Dendritic Tip Growth Approximation
                #define a3 0.00000149           // Estimated Parameter for Dendritic Tip Growth Approximation

                #define liquiddensity 2370.0
                #define soliddensity 2535.0
                #define viscosity 0.001

                #define Cp 921.5        // J/kg/K
                #define Lf 372000.0      // J/kg

                #define ml -7.7        // Liquidus Slope
                #define kpart 0.117      // Partition Coefficient
                #define Tm 933.0        // Melt Temperature

                #define GibbsThompson 9.8e-8   // Gibb-Thomspon Coefficient
                #define DiffLiq 3.0e-9
                //double Kl(){return 176.0;}
                //double Ks(){return 200.0;}
                
        #endif



        #ifdef Co7

                #define Co 7.0

                #define a2 0.00000290           // Estimated Parameter for Dendritic Tip Growth Approximation
                #define a3 0.00000149           // Estimated Parameter for Dendritic Tip Growth Approximation
                #define Tliq 891.0         // Liquidus Temperature for the given composition (7% w Si) (Kelvin)
                #define Tf 936.5
                #define Teut 850.0

                #define liquiddensity 2370.0
                #define soliddensity 2452.0
                //#define soliddensity 2000.0
                #define viscosity 0.001

                #define Cp 1126.0
                #define Lf 387400.0      // J/kg

                #define ml -6.5         // Liquidus Slope
                #define Tm 933.0        // Melt Temperature
                #define kpart 0.13      // Partition Coefficient

                #define GibbsThompson 1.96e-7   // Gibb-Thomspon Coefficient
                #define DiffLiq 6.45e-9      // Solute Diffusivity




        #endif




#endif

#ifdef AlCu

        #ifdef Co3

                #define Co 3.0

                #define Tliq 922.89         // Liquidus Temperature for the given composition (7% w Si) (Kelvin)
                #define Tf 933.0
                #define Teut 821.0         // Eutetic Temperature (Estimated Value) (ºC)

                #define a2 0.00000290           // Estimated Parameter for Dendritic Tip Growth Approximation
                #define a3 0.00000149           // Estimated Parameter for Dendritic Tip Growth Approximation

                #define liquiddensity 2370.0
                #define soliddensity 2633.2
                #define viscosity 0.001

                #define Cp 1291.2        // J/kg/K
                #define Lf 372000.0      // J/kg (Lf constant hypothesis)

                #define ml -3.37        // Liquidus Slope
                #define kpart 0.17      // Partition Coefficient
                #define Tm 933.0        // Melt Temperature

                #define GibbsThompson 9.8e-8   // Gibb-Thomspon Coefficient
                #define DiffLiq 5.0e-9
                //double Kl(){return 176.0;}
                //double Ks(){return 200.0;}

        #endif

        #ifdef Co1

                #define Co 1.0

                #define Tliq 933.0        // Liquidus Temperature for the given composition (7% w Si) (Kelvin)
                #define Tf 936.37
                #define Teut 821.0         // Eutetic Temperature (Estimated Value) (ºC)

                #define a2 0.00000290           // Estimated Parameter for Dendritic Tip Growth Approximation
                #define a3 0.00000149           // Estimated Parameter for Dendritic Tip Growth Approximation

                #define liquiddensity 2370.0
                #define soliddensity 2633.2
                #define viscosity 0.001

                #define Cp 1291.2        // J/kg/K
                #define Lf 372000.0      // J/kg (Lf constant hypothesis)

                #define ml -3.37        // Liquidus Slope
                #define kpart 0.14      // Partition Coefficient
                #define Tm 933.0        // Melt Temperature

                #define GibbsThompson 9.8e-8   // Gibb-Thomspon Coefficient
                #define DiffLiq 5.0e-9
                //double Kl(){return 176.0;}
                //double Ks(){return 200.0;}

        #endif

        #ifdef Co10

                #define Co 10.0

                #define Tliq 908.0         // Liquidus Temperature for the given composition (7% w Si) (Kelvin)
                #define Tf 941.7
                #define Teut 821.0         // Eutetic Temperature (Estimated Value) (ºC)

                #define a2 0.00000290           // Estimated Parameter for Dendritic Tip Growth Approximation
                #define a3 0.00000149           // Estimated Parameter for Dendritic Tip Growth Approximation

                #define liquiddensity 2582.0
                #define soliddensity 2589.0
                #define viscosity 0.001

                #define Cp 1291.2        // J/kg/K
                #define Lf 372000.0      // J/kg (Lf constant hypothesis)

                #define ml -3.37        // Liquidus Slope
                #define kpart 0.14      // Partition Coefficient
                #define Tm 933.0        // Melt Temperature

                #define GibbsThompson 9.8e-8   // Gibb-Thomspon Coefficient
                #define DiffLiq 5.0e-9
                //double Kl(){return 176.0;}
                //double Ks(){return 200.0;}

        #endif

#endif

#endif
