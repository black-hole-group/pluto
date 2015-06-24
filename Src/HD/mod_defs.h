/* ****************************************************************
     Set label, indexes and basic prototyping for the Hydro (HD)
     module.

     We make extra vector components, when not needed, point 
     to the last element (255) of the array stored by startup.c.  
   **************************************************************** */

enum {
 #if COMPONENTS == 1
  RHO, MX1, 
  #if HAVE_ENERGY
   ENG, PRS = ENG,
  #endif
  MX2 = 255, MX3 = 255,
 #elif COMPONENTS == 2
  RHO, MX1, MX2, 
  #if HAVE_ENERGY
   ENG, PRS = ENG,
  #endif
  MX3 = 255,
 #elif COMPONENTS == 3
  RHO, MX1, MX2, MX3,
  #if HAVE_ENERGY
   ENG, PRS = ENG,
  #endif
 #endif
 VX1 = MX1, VX2 = MX2, VX3 = MX3, 
 
/* -- backward compatibility -- */

 DN = RHO, 
 #if HAVE_ENERGY
  PR = PRS, EN = ENG,
 #endif
 VX = VX1, VY = VX2, VZ = VX3,
 MX = MX1, MY = MX2, MZ = MX3
};

#define NFLX (COMPONENTS + (HAVE_ENERGY == YES ? 2:1))

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

#endif

#if GEOMETRY == POLAR 

 #define iVR    VX1
 #define iVPHI  VX2
 #define iVZ    VX3

 #define iMR    MX1
 #define iMPHI  MX2
 #define iMZ    MX3

#endif

#if GEOMETRY == SPHERICAL 

 #define iVR    VX1
 #define iVTH   VX2
 #define iVPHI  VX3

 #define iMR    MX1
 #define iMTH   MX2
 #define iMPHI  MX3

#endif

/* *************************************************
     Label the different waves in increasing order 
     following the number of vector components.
   ************************************************* */

enum KWAVES {
 KSOUNDM, KSOUNDP
 #if HAVE_ENERGY
  , KENTRP
 #endif
};

/* ***********************************************************
                   Prototyping goes here          
   *********************************************************** */

int  ConsToPrim   (double **, double **, int, int, unsigned char *);
void Eigenvalues (double **, double *, double **, int, int);
void PrimEigenvectors (double *, double, double, double *, double **, double **);
void ConsEigenvectors (double *, double *, double, 
                       double **, double **, double *);

void Flux      (double **, double **, double *, double **, double *, int, int);
void HLL_Speed (double **, double **, double *, double *, 
                double *, double *, int, int);
void MaxSignalSpeed (double **, double *, double *, double *, int, int);
void PrimToCons   (double **, double **, int, int);
void PrimRHS    (double *, double *, double, double, double *);
void PrimSource (const State_1D *, int, int, 
                 double *, double *, double **, Grid *);

Riemann_Solver TwoShock_Solver, LF_Solver, Roe_Solver, HLL_Solver, 
               HLLC_Solver, RusanovDW_Solver;
Riemann_Solver AUSMp_Solver;


