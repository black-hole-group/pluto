/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief PLUTO main header file.

  Contains basic macro definitions, structure definitions and global 
  variable declarations used by the code.

  \author A. Mignone (mignone@ph.unito.it)
  \date   Aug 16, 2012
*/
/* ///////////////////////////////////////////////////////////////////// */
#ifndef PLUTO_H
#define PLUTO_H

#define PLUTO_VERSION  "4.1"

#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

#define YES      1
#define NO       0
#define DEFAULT -1
#define TRUE    YES
#define FALSE    NO

/* ---- Geometry Labels ( > 0) ----  */

#define CARTESIAN    1
#define CYLINDRICAL  2
#define POLAR        3
#define SPHERICAL    4

#define UNIFORM_GRID             1
#define STRETCHED_GRID           2
#define LOGARITHMIC_INC_GRID     3
#define LOGARITHMIC_DEC_GRID     4

/* ---- Equation of state (EoS) labels  ----  */

#define IDEAL         1
#define PVTE_LAW      2
#define TAUB          3
#define BAROTROPIC    4
#define ISOTHERMAL    5

/* ---- Time stepping labels ----  */

#define EULER                      1
#define HANCOCK                    2
#define CHARACTERISTIC_TRACING     3
#define RK2                        5
#define RK3                        6
#define RK_MIDPOINT                7

#define EXPLICIT             1 /* -- just a number different from 0 !!!  -- */
#define SUPER_TIME_STEPPING  2 /* -- just a number different from EXPLICIT -- */ 
#define RK_CHEBYSHEV         4  

/* ---- Operator step labels ---- */

#define HYPERBOLIC_STEP  1  
#define PARABOLIC_STEP   2
#define SOURCE_STEP      3

/* ----   Output labels  ---- */

#define DBL_OUTPUT      1
#define FLT_OUTPUT      2
#define VTK_OUTPUT      3
#define DBL_H5_OUTPUT   4
#define FLT_H5_OUTPUT   5
#define TAB_OUTPUT      6
#define PPM_OUTPUT      7
#define PNG_OUTPUT      8

#define VTK_VECTOR  5  /* -- any number but NOT 1  -- */

/*! The maximum number of output file formats is fixed to 11 so that the 
    size of runtime structure (defined below) is 64 bytes. 
    This should prevent, for some compilers, attempts to change the 
    alignment of the structure and therefore troubleshooting when restarting 
    from files written on different architectures.                    */
#define MAX_OUTPUT_TYPES 11 
#define MAX_OUTPUT_VARS  128

/* ----  Cooling labels ----  */

#define POWER_LAW    3
#define MINEq        4
#define SNEq         5
#define TABULATED    6
#define H2_COOL      7

/* ---- Physics modules labels ----  */

#define ADVECTION  1
#define HD         2
#define RHD        3
#define MHD        4
#define RMHD       5

    /*  ----  SET LABELS FOR DIV.B REMOVAL  ----  
        If you move them to the MHD header, 
        definitions.h (which is included before)
        cannot correctly use them                */
        
#define NONE                   0
#define EIGHT_WAVES            1
#define DIV_CLEANING           2
#define CONSTRAINED_TRANSPORT  3

   /*  ----  SET LABELS FOR BODY_FORCE  ----
      Please do not change them since they are
      used in bitwise operations                */
   
#define VECTOR     4   /* corresponds to  100 in binary  */
#define POTENTIAL  8   /* corresponds to 1000 in binary  */

/* ---- Boundary condition labels  ----  */

#define OUTFLOW          1  /* any number except 0 !! */
#define REFLECTIVE       2 
#define AXISYMMETRIC     3
#define EQTSYMMETRIC     4
#define PERIODIC         5
#define SHEARING         6
#define USERDEF          7

#define X1_BEG           101
#define X1_END           102
#define X2_BEG           103
#define X2_END           104
#define X3_BEG           105
#define X3_END           106

/* ---- LABELS FOR IMAGE SLICING ---- */

#define X12_PLANE       3  
#define X13_PLANE       5
#define X23_PLANE       6

/*! \name Bit flag labels.
    The following macros define the bits that can be turned on or off
    in an unsigned char (= 1 byte = 8 bits) variable. 
    Different bit flags allow to enable or disable certain actions in 
    a given cell at different points in the code, see also flag.c.
    The 3D unsigned char \c ***flag array is used for bookeeping, in each zone 
    (i,j,k), which bits are actually switched on or off.
    A simple bitwise operation is used to enable a flag, e.g., 
    <tt> flag[k][j][i] |= FLAG_XXX </tt>.
    For instance, by turning the ::FLAG_HLL bit on, we have
    <tt> flag = 00000100 </tt>, while by also enabling the ::FLAG_SPLIT_CELL
    one has <tt> flag = 00010100 </tt> and so on.
    Individual bits can also be turned off with the complement operator

    Individual bits can also be turned off, e.g., <tt> flag &= ~FLAG_XXX; </tt>

*/
/**@{ */
#define FLAG_MINMOD      1  /**< Reconstruct using MINMOD limiter. */
#define FLAG_FLAT        2  /**< Reconstruct using FLAT limiter.   */
#define FLAG_HLL         4  /**< Use HLL Riemann solver. */
#define FLAG_ENTROPY     8  /**< Update pressure using entropy equation. */
#define FLAG_SPLIT_CELL  16 /**< Zone is covered by a finer level (AMR only). */
#define FLAG_INTERNAL_BOUNDARY   32  /**< Zone belongs to an internal boundary
                                          region and should be excluded from 
                                          being updated in time              */
#define FLAG_BIT7        64
#define FLAG_BIT8       128  
/**@} */

#define PRS_FAIL  1
#define ENG_FAIL  2
#define RHO_FAIL  4

#define IDIR     0     /*   This sequence (0,1,2) should */
#define JDIR     1     /*   never be changed             */
#define KDIR     2     /*                                */
#define ALL_DIR -1

/* -- location of a variable inside the cell -- */

#define CENTER  0
#define X1FACE  1
#define X2FACE  2
#define X3FACE  3

#define CELL_CENTER    50  /* really needed ? */
#define FACE_CENTER    51
#define EDGE_CENTER    52

        /*  ----  SET LABELS FOR INTERPOLATION  ----   */

#define FLAT              1
#define LINEAR            2
#define CENO3             3
#define PARABOLIC         4
#define LINEAR_MULTID     5
#define MP5               6
#define LimO3             7
#define WENO3             8

#define WENO3_FD             103
#define WENO5_FD             105
#define WENOZ_FD             106
#define WENO7_FD             107
#define MP5_FD               125
#define PPM_FD               140
#define LIMO3_FD             300

#define ONED   1
#define MULTID 3

/* ---- limiter labels ---- */

#define FLAT_LIM          1
#define MINMOD_LIM        2 
#define VANALBADA_LIM     3
#define OSPRE_LIM         4
#define UMIST_LIM         5
#define VANLEER_LIM       6
#define MC_LIM            7
#define FOURTH_ORDER_LIM  8

/*! \name Physical constants in c.g.s units.
     The following set of macros express some useful physical constants
     in c.g.s units (erg, cm and sec). Values have been taken from
     http://physic.nist.gov/cuu/Constants/index.html
*/
/**@{ */
#define CONST_AH      1.008              /**< Atomic weight of Hydrogen  */
#define CONST_AHe     4.004              /**< Atomic weight of Helium  */
#define CONST_AZ      30.0               /**< Mean atomic weight of heavy elements */
#define CONST_amu     1.66053886e-24     /**<  Atomic mass unit.          */
#define CONST_au      1.49597892e13      /**<  Astronomical unit.         */
#define CONST_c       2.99792458e10      /**<  Speed of Light.            */
#define CONST_eV      1.602176463158e-12 /**<  Electron Volt in erg.      */
#define CONST_G       6.6726e-8          /**<  Gravitational Constant.    */
#define CONST_h       6.62606876e-27     /**<  Planck Constant.           */
#define CONST_kB      1.3806505e-16      /**<  Boltzmann constant.        */
#define CONST_ly      0.9461e18          /**<  Light year.                */
#define CONST_mp      1.67262171e-24     /**<  Proton mass.               */
#define CONST_mn      1.67492728e-24     /**<  Neutron mass.              */
#define CONST_me      9.1093826e-28      /**<  Electron mass.             */
#define CONST_mH      1.6733e-24         /**<  Hydrogen atom mass.        */
#define CONST_Msun    2.e33              /**<  Solar Mass.                */
#define CONST_Mearth  5.9736e27          /**<  Earth Mass.                */
#define CONST_NA      6.0221367e23       /**<  Avogadro Contant.          */
#define CONST_pc      3.0856775807e18    /**<  Parsec.                    */
#define CONST_PI      3.14159265358979   /**<  \f$ \pi \f$.               */
#define CONST_Rearth  6.378136e8         /**<  Earth Radius.              */
#define CONST_Rsun    6.96e10            /**<  Solar Radius.              */
#define CONST_sigma   5.67051e-5         /**<  Stephan Boltmann constant. */
#define CONST_sigmaT  6.6524e-25         /**<  Thomson Cross section.    */
/**@} */

/* **********************************************************
    Start including header files here
   ********************************************************** */

#include "definitions.h"   /* Problem-dependent header file  */

/* *********************************************************************
    Define unit density, length and velocity if not done before
   ********************************************************************* */

#ifndef UNIT_DENSITY
 #define UNIT_DENSITY (CONST_mp)  /**< Unit density in gr/cm^3. */
#endif

#ifndef UNIT_LENGTH
 #define UNIT_LENGTH   (CONST_au)  /**< Unit Legnth in cm. */
#endif

#ifndef UNIT_VELOCITY
 #define UNIT_VELOCITY (1.e5)  /**< Unit velocity in cm/sec. */
#endif

/*! Define the conversion constant between dimensionless 
    temperature prs/rho and physical temperature in Kelvin,
    T = (prs/rho)*KELVIN*mu                                   */
#define KELVIN (UNIT_VELOCITY*UNIT_VELOCITY*CONST_amu/CONST_kB) 

/* *********************************************************************
    set to NO all undefined identifier that may appear in the code but
    have not been defined in definitions.h
   ********************************************************************* */

#ifndef RESISTIVE_MHD 
 #define RESISTIVE_MHD   NO
#endif

#ifndef THERMAL_CONDUCTION
 #define THERMAL_CONDUCTION NO
#endif

#ifndef VISCOSITY
 #define VISCOSITY NO
#endif

#ifndef INCLUDE_PARTICLES
 #define INCLUDE_PARTICLES NO
#endif

#ifndef UPDATE_VECTOR_POTENTIAL
 #define UPDATE_VECTOR_POTENTIAL  NO
#endif

#ifndef BACKGROUND_FIELD
 #define BACKGROUND_FIELD NO
#endif

#ifndef USE_FOUR_VELOCITY
 #define USE_FOUR_VELOCITY NO
#endif

#ifndef CHAR_LIMITING
 #define CHAR_LIMITING  NO
#endif

#ifndef ARTIFICIAL_VISCOSITY 
 #define ARTIFICIAL_VISCOSITY  NO
#endif

#ifndef ROTATING_FRAME
 #define ROTATING_FRAME NO
#endif

#ifdef CH_SPACEDIM
 #define CHOMBO  1

 #ifndef CHOMBO_LOGR
  #define CHOMBO_LOGR NO
 #endif

 #ifndef CHOMBO_CONS_AM
  #if (GEOMETRY == CYLINDRICAL) && (COMPONENTS == 3) && (CHOMBO_EN_SWITCH == YES)
   #define CHOMBO_CONS_AM YES
  #elif (GEOMETRY == SPHERICAL) && (COMPONENTS == 3) && (CHOMBO_EN_SWITCH == YES)
   #define CHOMBO_CONS_AM YES
  #elif (GEOMETRY == POLAR) && (COMPONENTS > 1) && (CHOMBO_EN_SWITCH == YES)
   #define CHOMBO_CONS_AM YES
  #else
   #define CHOMBO_CONS_AM NO
  #endif
 #endif

 #if CHOMBO_CONS_AM == YES
  #define CHOMBO_NDV  2
 #else
  #define CHOMBO_NDV  1
 #endif

#endif

/* ******************************************************
    Debug switches
   ****************************************************** */
 
/* -- CHECK_EIGENVECTORS: used in eigenv.c in HD/, MHD/, RHD/
      to check orthogonality and the correctness through 
      the relation the A = L*\Lambda*R  -- */
#define CHECK_EIGENVECTORS      NO

/* -- CHECK_CONSERVATIVE_VAR: used in RHD/mappers.c to 
      check that conservative vars are physical -- */

#define CHECK_CONSERVATIVE_VAR  NO

/* -- CHECK_DIVB_CONDITION: used in MHD/CT/ct.c 
      to check if div.B = 0 -- */

#ifndef CHECK_DIVB_CONDITION
 #define CHECK_DIVB_CONDITION   NO
#endif

/* -- Shortcut for CTU -- */

#if (TIME_STEPPING == HANCOCK) || (TIME_STEPPING == CHARACTERISTIC_TRACING)
 #if DIMENSIONAL_SPLITTING == NO
  #define CTU      1    /* -- Corner Transport Upwind method of Colella -- */
 #endif
#endif

/* ------------------------------------------------------------
    the GET_MAX_DT switch determines how the time step is
    computed. For pure advection, setting it to YES will force
    PLUTO to compute dt in the old way, i.e., by taking the
    maximum.
    When set to NO, the time step is computed only during the
    predictor step and, for UNSPLIT RK schemes, it will be
    calculated as the average over dimensions resulting in 
    slightly larger time increments.
   ------------------------------------------------------------ */

#if ((TIME_STEPPING == RK2) || (TIME_STEPPING == RK3) || \
     (TIME_STEPPING == RK_MIDPOINT)) && DIMENSIONAL_SPLITTING == NO
 #define GET_MAX_DT    NO
#else
 #define GET_MAX_DT    YES
#endif

/* *********************************************************************
    Diffusion operators: PARABOLIC_FLUX is the bitwise OR
    combination of all operators, each being either one of 
    NO, EXPLICIT (1st bit), STS (2nd bit). 
    It can take the following values

      00   --> no diffusion operator is being used
      01   --> there's at least one explicit diffusion operator and
               no sts.
               
      10   --> there's at least one sts diffusion operator and
               no explicit one.
      11   --> mixed: there is at least one explicit and sts operator
   ********************************************************************* */

#define PARABOLIC_FLUX (RESISTIVE_MHD|THERMAL_CONDUCTION|VISCOSITY)

/* ********************************************************
    Include more header files
   ******************************************************** */

#ifdef PARALLEL      /* Only for parallel computations on static grid */
 #include <al.h>
#endif
#ifdef CH_MPI       /* Include mpi.h for parallel Chombo, in order to */
 #include <mpi.h>   /* use MPI_Abort function in the QUIT_PLUTO macro */
#endif
#include "macros.h"  /* Function-like macro header file */
#include "structs.h" /* Structure declaration header file */

/* *****************************************************
     Recurrent types
     Note: when using Finite Difference Schemes, the
          "Riemann Solver" function computes the fluxes
           with high order interpolants. 
   ***************************************************** */

typedef double real;
typedef void Riemann_Solver (const State_1D *, int, int, double *, Grid *);
typedef void Limiter        (double *, double *, double *, int, int, Grid *);
typedef double Reconstruct  (double *, double, int);
typedef double ****Data_Arr;

#ifndef EOS
 #define EOS  -1
#endif

/* -- Set HAVE_ENERGY to YES if energy equation exists -- */
#if (EOS == IDEAL) || (EOS == PVTE_LAW) || (EOS == TAUB)  
  #define HAVE_ENERGY  YES
#else
  #define HAVE_ENERGY  NO
#endif
#define HAVE_PRESSURE  HAVE_ENERGY

/* ********************************************************
    Include yet more header files
   ******************************************************** */

#include "mod_defs.h"  /* Include physics header file (search path is set
                          in the makefile) */
#ifdef SHEARINGBOX
 #include "MHD/ShearingBox/shearingbox.h"
#endif

#if COOLING != NO    
 #include "cooling.h"               /* Include cooling header file */
#endif
#define RHOE PRS

#if THERMAL_CONDUCTION != NO 
 #include "Thermal_Conduction/tc.h" /* Thermal conduction header file */
#endif

#if VISCOSITY != NO
 #include "Viscosity/viscosity.h"   /* Viscosity header file */
#endif

#ifdef FARGO
 #include "Fargo/fargo.h"           /* FARGO header file */
#endif

#include "States/plm_coeffs.h"      /* PLM header file */
#if INTERPOLATION == PARABOLIC 
 #include "States/ppm_coeffs.h"     /* PPM header file */
#endif

#include "Math_Tools/math_tools.h"  /* Math tools header file */

/* *********************************************************************
    Define total number of variables to be integrated;
    this includes:

     NFLX     = Number of equations defining the system of 
                conservation laws. For example, for the 
                HD module, it consists of density, momentum and energy.
                It is defined in the physics module header file mod_defs.h.

     NTRACER  = Number of user-defined tracers; defined in the problem
                directory header file definitions.h

     NIONS = Number of chemical species; defined in the cooling modules 
             cooling.h, if present.

  nv = 0...NFLX - 1; NFLX...NFLX+NIONS-1; TRC...TRC + NTRACER-1; ENTR; 
       <----------> <------------------> <-------------------->  
          NFLX              NIONS               NTRACER              

                    <----------------------------------------------->
                                         NSCL
                                         
       <------------------------------------------------------------>
                                   NVAR

   ********************************************************************* */

#ifndef ENTROPY_SWITCH
 #define ENTROPY_SWITCH  NO 
#endif 

#ifndef NIONS
 #define NIONS 0
#endif

#define TRC  (NFLX + NIONS)
#define TR    TRC    /* for backward compatibility */
#define NSCL (NTRACER + NIONS + ENTROPY_SWITCH)

#if ENTROPY_SWITCH == YES
 #define ENTR  (TRC + NTRACER)
#else
 #if HAVE_ENERGY
  #define ENTR (ENG)
 #endif
#endif

/* -- define total number of variables to be evolved in time -- */

#define NVAR (NFLX + NSCL)

/* ----------------------------------------------------------
    Include module header files: EOS
    [This section should be placed before, but NVAR 
     wouldn't be defined. Need to fix this at some point]
   ---------------------------------------------------------- */

#include "eos.h"
#include "prototypes.h"

/* *****************************************************
     Declare global variables
   ***************************************************** */

 extern int SZ;
 extern int SZ_stagx;
 extern int SZ_stagy;
 extern int SZ_stagz;
 extern int SZ_char;
 extern int SZ_float;
 extern int SZ_Float_Vect;
 extern int SZ_rgb;
 extern int SZ_short;
 extern int prank;

extern long int IBEG, IEND, JBEG, JEND, KBEG, KEND;
extern long int NX1, NX2, NX3;
extern long int NX1_TOT, NX2_TOT, NX3_TOT;
extern long int NMAX_POINT;

extern int VXn, VXt, VXb;
extern int MXn, MXt, MXb;
extern int BXn, BXt, BXb;

extern int g_i, g_j, g_k;

extern int g_dir;
extern int g_maxRiemannIter;
extern int g_maxRootIter;
extern long int g_usedMemory;
extern long int g_stepNumber;
extern int      g_intStage;
extern int      g_operatorStep;

extern double g_maxCoolingRate, g_minCoolingTemp;

extern double g_smallDensity, g_smallPressure;

extern double g_time, g_dt;
extern double g_maxMach;
#if ROTATING_FRAME
 extern double g_OmegaZ;
#endif

extern double g_domBeg[3], g_domEnd[3];

extern double g_inputParam[32];
#if EOS == IDEAL
 extern double g_gamma;
#elif EOS == ISOTHERMAL
 extern double g_isoSoundSpeed;
#endif

#ifdef CHOMBO
 extern double glm_ch_max, glm_ch_max_loc, g_coeff_dl_min;
 extern double g_level_dx;
 extern double g_x2stretch, g_x3stretch;
 extern int    glm_is_defined;
 #if GEOMETRY == CARTESIAN
  extern double g_stretch_fact;
 #endif
#endif

/* ---- Maximum grid size for allocating static arrays ---- */

#ifdef CHOMBO

 #define NX1_MAX   NMAX_POINT
 #if DIMENSIONS == 1
  #define NX2_MAX   1
  #define NX3_MAX   1
 #elif DIMENSIONS == 2
  #define NX2_MAX   NMAX_POINT
  #define NX3_MAX   1
 #else
  #define NX2_MAX   NMAX_POINT
  #define NX3_MAX   NMAX_POINT
 #endif

#else

 #define NX1_MAX   NX1_TOT
 #define NX2_MAX   NX2_TOT
 #define NX3_MAX   NX3_TOT

#endif

#endif /* PLUTO_H */
