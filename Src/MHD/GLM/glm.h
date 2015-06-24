/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief  Header file for GLM Divergence Cleaning

  Contains function prototypes and global variable declaration
  for the GLM formulation to control the divergence-free condition 
  of magnetic field.

  \b References
     - "A Second-order unsplit Godunov scheme for cell-centered MHD:
        The CTU-GLM scheme"\n
        Mignone \& Tzeferacos, JCP (2010) 229, 2117

     - "High-order conservative finite difference GLM-MHD scheme for
        cell-centered MHD"\n
        Mignone, Tzeferacos \& Bodo, JCP (2010) 229, 5896
        
  \authors A. Mignone (mignone@ph.unito.it)\n
           P. Tzeferacos (petros.tzeferacos@ph.unito.it)
  \date    Sept 23, 2014
*/
/* ///////////////////////////////////////////////////////////////////// */
#define GLM_MHD

#ifndef EGLM
 #define EGLM   NO  /**< The EGLM macro may be turned to YES to enable
                         the extended GLM formalism. 
                         Although it breaks conservation of momentum and
                         energy, it has proven to be more robust in treating
                         low-beta plasma. */
#endif

#define COMPUTE_DIVB  NO

/* with chombo, COMPUTE_DIVB must be 
   disabled or a segfault will occur */

#ifdef CHOMBO
 #undef COMPUTE_DIVB  
 #define COMPUTE_DIVB NO
#endif

extern double glm_ch; /**< The propagation speed of divergence error. */
    
void  GLM_Solve (const State_1D *, double **, double **, int, int, Grid *);
void  GLM_SolveNEW (const State_1D *state, int beg, int end, Grid *grid);
void  GLM_Init      (const Data *, const Time_Step *, Grid *);
void  GLM_Source (const Data_Arr, double, Grid *);
void  EGLM_Source (const State_1D *, double, int, int, Grid *);

#if COMPUTE_DIVB == YES
 void GLM_ComputeDivB(const State_1D *state, Grid *grid);
 double ***GLM_GetDivB(void);
#endif

