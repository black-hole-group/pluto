/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Set labels, indexes and prototypes for the MHD module.

  Contains basic macro definitions, structure definitions and global
  variable declarations used by the MHD module.

  \author A. Mignone (mignone@ph.unito.it)
  \date   Sep, 10, 2012
*/
/* ///////////////////////////////////////////////////////////////////// */


/* ****************************************************************
     We make extra vector components, when not needed, point 
     to the last element (255) of the array stored by startup.c.  
   **************************************************************** */

enum {

 #if COMPONENTS == 1

  RHO, MX1, BX1, 
  #if HAVE_ENERGY
   ENG, PRS = ENG,
  #endif
  #if MHD_FORMULATION == DIV_CLEANING
   PSI_GLM,
  #endif
  MX2 = 255, BX2 = 255, MX3 = 255, BX3 = 255,

 #elif COMPONENTS == 2

  RHO, MX1, MX2, BX1, BX2,
  #if HAVE_ENERGY
   ENG, PRS = ENG,
  #endif
  #if MHD_FORMULATION == DIV_CLEANING
   PSI_GLM,
  #endif
  MX3 = 255, BX3 = 255,

 #elif COMPONENTS == 3

  RHO, MX1, MX2, MX3, BX1, BX2, BX3, 
  #if HAVE_ENERGY
   ENG, PRS = ENG,
  #endif
  #if MHD_FORMULATION == DIV_CLEANING
   PSI_GLM,
  #endif

 #endif

 VX1 = MX1, VX2 = MX2, VX3 = MX3, 

/* -- backward compatibility -- */

 DN = RHO,
 #if HAVE_ENERGY
  PR = PRS, EN = ENG,
 #endif
 VX = VX1, VY = VX2, VZ = VX3,
 MX = MX1, MY = MX2, MZ = MX3,
 BX = BX1, BY = BX2, BZ = BX3

};

#define NFLX (2*COMPONENTS + (HAVE_ENERGY == YES ? 2:1) \
                           + (MHD_FORMULATION == DIV_CLEANING))

/* ********************************************************************* */
/*! Label the different waves in increasing order 
    following the number of vector components.

    \b IMPORTANT: the KPSI_GLMM & KPSI_GLMP modes are 
                 present only in the MHD-GLM formulation.
                 We keep them at the END of the enumeration
                 so we can skip them in unnecessary loops.
                 Please do NOT change them !
   ********************************************************************* */

enum KWAVES {
 KFASTM, KFASTP
 #if HAVE_ENERGY
  , KENTRP
 #endif

 #if MHD_FORMULATION != DIV_CLEANING
  , KDIVB
 #endif

 #if COMPONENTS >= 2
  , KSLOWM, KSLOWP
  #if COMPONENTS == 3
   , KALFVM, KALFVP
  #endif
 #endif

 #if MHD_FORMULATION == DIV_CLEANING  
  , KPSI_GLMM, KPSI_GLMP 
 #endif
};

/*! \name Vector Potential Labels 
    These may only be used in the STARTUP / INIT  functions.
    They're convenient in obtaining a discretization that preserve 
    the divergence-free condition (for staggered field) or if you simply
    wish to initialize the magnetic field from the vector potential.     */
/**@{ */
#define   AX1  (NVAR + 1)
#define   AX2  (NVAR + 2)
#define   AX3  (NVAR + 3)
/**@} */

#define AX  AX1  
#define AY  AX2
#define AZ  AX3

/* *************************************************
     Now define more convenient and user-friendly 
     pointer labels for geometry setting      
   ************************************************* */

#if GEOMETRY == CYLINDRICAL 

 #define iVR    VX1
 #define iVZ    VX2
 #define iVPHI  VX3

 #define iMR    MX1
 #define iMZ    MX2
 #define iMPHI  MX3

 #define iBR    BX1
 #define iBZ    BX2
 #define iBPHI  BX3

#endif

#if GEOMETRY == POLAR 

 #define iVR    VX1
 #define iVPHI  VX2
 #define iVZ    VX3

 #define iMR    MX1
 #define iMPHI  MX2
 #define iMZ    MX3

 #define iBR    BX1
 #define iBPHI  BX2
 #define iBZ    BX3

#endif

#if GEOMETRY == SPHERICAL 

 #define iVR     VX1
 #define iVTH    VX2
 #define iVPHI   VX3

 #define iMR    MX1
 #define iMTH   MX2
 #define iMPHI  MX3

 #define iBR    BX1
 #define iBTH   BX2
 #define iBPHI  BX3

#endif

/* ***********************************************************
    \cond REPEAT_FUNCTION_DOCUMENTATION_IN_HEADER_FILES 
    Function prototyping
   *********************************************************** */

void BackgroundField (double x1, double x2, double x3, double *B0);

#if BACKGROUND_FIELD == YES
 double **GetBackgroundField (int, int, int, Grid *);
#endif 

int  ConsToPrim   (double **, real **, int , int, unsigned char *);
void Eigenvalues (double **, double *, double **, int, int);

void PrimEigenvectors (double *, double, double, double *, double **, double **);
void ConsEigenvectors (double *, double *, double,
                       double **, double **, double *);

void Flux      (double **, double **, double *, double **, double **,
                double *, int, int);
void HLL_Speed (double **, double **, double *, double *, double **, 
                double *, double *, int, int);
void MaxSignalSpeed (double **, double *, double *, double *, double **, int, int);
void PrimToCons   (double **, double **, int, int);
void PrimRHS    (double *, double *, double, double, double *);
void PrimSource (const State_1D *, int, int, 
                 double *, double *, double **, Grid *);

#if MHD_FORMULATION == EIGHT_WAVES

 void Roe_DivBSource (const State_1D *, int, int, Grid *);
 void HLL_DivBSource (const State_1D *, double **, int, int, Grid *);

#elif MHD_FORMULATION == DIV_CLEANING

 #include "MHD/GLM/glm.h"

#elif MHD_FORMULATION == CONSTRAINED_TRANSPORT

 #include "MHD/CT/ct.h"

#endif

Riemann_Solver HLL_Solver, HLLC_Solver, HLLD_Solver;
Riemann_Solver LF_Solver, Roe_Solver, MUSTA_Solver;

#if RESISTIVE_MHD != NO
 #include "Resistive/res.h"
#endif

#ifdef SHEARINGBOX
 #include "MHD/ShearingBox/shearingbox.h"
#endif
/* \endcond */
