#include "pluto.h"

/* ********************************************************************* */
void Init (double *us, double x1, double x2, double x3)
/*
 *  Taylor-Couette Flow. 
 *
 *  Relevant reading:
 *
 *  Taylor G.I., 
 *  "Stability of a viscous liquid contained between two rotating cylinders"
 *  1923, Philos. Trans. R. Soc. London, Ser. A, 223, 289
 * 
 *  Recktenwald A., L\"ucke M., M\"uller H.W.,  
 *  "Taylor vortex formation in axial through-flow: Linear and 
 *  weakly nonlinear analysis" 
 *  1993, Phys. Rev. E, 48, 4444 
 *  
 *  The fluid is rotating between two cylinders situated at
 *  R_{INT} = 0.9 (rotating at \Omega = 1.0) and R_{EXT} = 1.0 (\Omega = 0). 
 *  Viscous effects are controlled by a Reynolds number defined in 
 *  eta_visc.c as Re = \Omega * R_{INT} * (R_{EXT}-R_{INT}) * \rho /eta1_visc
 *  For Reynolds numbers Re > Re_{critical}~130  vortices are formed, the axial distribution 
 *  of which is controlled by the wave-number of the initial perturbation 
 *  \kappa = 2 \pi / \lambda.
 *  For small Reynolds numbers, viscosity suppresses the vortex formation.
 *  (See included .jpeg files)
 *
 *  Written by: Petros Tzeferacos (petros.tzeferacos@ph.unito.it)
 *  Last modified: 04 June 2010
 *
 *
 *********************************************************************** */
{
  double eta = g_inputParam[R_INT]/g_inputParam[R_EXT];
  double lambda = (g_inputParam[R_EXT]-g_inputParam[R_INT]); 
  double A = -g_inputParam[OMEGA]*eta*eta/(1.-eta*eta);
  double B =  g_inputParam[OMEGA]*g_inputParam[R_INT]*g_inputParam[R_INT]/(1.-eta*eta);
  double SPER = sin(2.*CONST_PI*x2/lambda);
  double CPER = cos(2.*CONST_PI*x2/lambda);
  
  
  
  us[RHO] = 1.0; 
  us[PRS] = 10.+ 0.5*(A*A*x1*x1 +4.*A*B*log(x1) - B*B/x1/x1) 
           + 0.01*0.5*(A*A*x1*x1 +4.*A*B*log(x1) - B*B/x1/x1)*CPER; 
  us[VX1] = 0.01*(A*x1 + B/x1)*CPER;
  us[VX2] = 0.01*(A*x1 + B/x1)*SPER;
  us[VX3] = A*x1 + B/x1 + 0.01*(A*x1 + B/x1)*CPER;
  us[TRC] = 1.0;

  #if PHYSICS == MHD || PHYSICS == RMHD

   us[BX1] = 0.0;
   us[BX2] = 0.0;
   us[BX3] = 0.0;

   #ifdef STAGGERED_MHD

    us[AX] = 0.0;
    us[AY] = 0.0;
    us[AZ] = 0.0;

   #endif

  #endif
}
/* ********************************************************************* */
void Analysis (const Data *d, Grid *grid)
/* 
 *
 *
 *********************************************************************** */
{

}
/* ********************************************************************* */
void UserDefBoundary (const Data *d, RBox *box, int side, Grid *grid) 
/*! 
 *  Assign user-defined boundary conditions.
 *
 * \param [in,out] d  pointer to the PLUTO data structure containing
 *                    cell-centered primitive quantities (d->Vc) and 
 *                    staggered magnetic fields (d->Vs, when used) to 
 *                    be filled.
 * \param [in] box    pointer to a RBox structure containing the lower
 *                    and upper indices of the ghost zone-centers/nodes
 *                    or edges at which data values should be assigned.
 * \param [in] side   specifies the boundary side where ghost zones need
 *                    to be filled. It can assume the following 
 *                    pre-definite values: X1_BEG, X1_END,
 *                                         X2_BEG, X2_END, 
 *                                         X3_BEG, X3_END.
 *                    The special value side == 0 is used to control
 *                    a region inside the computational domain.
 * \param [in] grid  pointer to an array of Grid structures.
 *
 *********************************************************************** */
{
  int   i, j, k, nv;
  real  *x1, *x2, *x3;

  x1 = grid[IDIR].x;
  x2 = grid[JDIR].x;
  x3 = grid[KDIR].x;

  if (side == 0) {    /* -- check solution inside domain -- */
    DOM_LOOP(k,j,i){};
  }

  if (side == X1_BEG){  /* -- X1_BEG boundary / Cylinder @ R = R_INT -- */
    X1_BEG_LOOP(k,j,i){
      d->Vc[RHO][k][j][i] =   d->Vc[RHO][k][j][2*IBEG - i - 1];
      d->Vc[PRS][k][j][i] =   d->Vc[PRS][k][j][2*IBEG - i - 1];
      d->Vc[VX1][k][j][i] = - d->Vc[VX1][k][j][2*IBEG - i - 1];
      d->Vc[VX2][k][j][i] =   d->Vc[VX2][k][j][2*IBEG - i - 1];
      d->Vc[VX3][k][j][i] =   g_inputParam[OMEGA]*x1[i];	    	    
    }
  }

  if (side == X1_END){  /* -- X1_END boundary / Cylinder @ R = R_EXT -- */
    X1_END_LOOP(k,j,i){
      d->Vc[RHO][k][j][i] =   d->Vc[RHO][k][j][2*IEND - i + 1];
      d->Vc[PRS][k][j][i] =   d->Vc[PRS][k][j][2*IEND - i + 1];
      d->Vc[VX1][k][j][i] = - d->Vc[VX1][k][j][2*IEND - i + 1];
      d->Vc[VX2][k][j][i] =   d->Vc[VX2][k][j][2*IEND - i + 1];
      d->Vc[VX3][k][j][i] =   0.0;	  	    
    }
  }
}

