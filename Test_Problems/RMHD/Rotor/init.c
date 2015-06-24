#include "pluto.h"

/* ********************************************************************* */
void Init (double *us, double x1, double x2, double x3)
/*
 *
 *
 *
 *********************************************************************** */
{
/* ---------------------------------------------------------------
     Set inital condition for the Rotor problem;
 
     Reference: 
 
     "On the divergence-free condition in Godunov-type schemes 
      for ideal MHD: the upwind constrained transport method"
  
     P. Londrillo, L. Del Zanna 
     JCP (2004), 195, 17
   --------------------------------------------------------------- */

  double r, r0, Bx, omega;

  g_gamma   = 5./3.;
  omega = 0.995;
  Bx    = 1.0;
  r0    = 0.1;

  #if GEOMETRY == CARTESIAN

   r = D_EXPAND(x1*x1, + x2*x2, + x3*x3);
   r = sqrt(r);

   us[RHO] = 1.0;
   us[VX1] = 0.0;
   us[VX2] = 0.0;
   us[PRS] = 1.0;
   us[BX1] = Bx;
   us[BX2] = 0.0;

   if (r <= r0) {
     us[RHO] = 10.0;
     us[VX1] = -omega*x2/r0;
     us[VX2] =  omega*x1/r0;
   }

   us[AX1] = us[AX2] = 0.0;
   us[AX3] = us[BX1]*x2;

  #elif GEOMETRY == POLAR

   r = x1;

   us[RHO] = 1.0;
   us[VX2] = 0.0;
   us[VX1] = 0.0;
   us[BX1] =  Bx*cos(x2);
   us[BX2] = -Bx*sin(x2);
   us[PRS] = 1.0;
  
   if (r <= r0) {
     us[RHO] = 10.0;
     us[VX2] = omega*r/r0;
   }

   us[AX1] = us[AX2] = 0.0;
   us[AX3] = - r*sin(x2)*Bx;

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
 * \param [in/out] d  pointer to the PLUTO data structure containing
 *                    cell-centered primitive quantities (d->Vc) and 
 *                    staggered magnetic fields (d->Vs, when used) to 
 *                    be filled.
 * \param [in] box    pointer to a RBox structure containing the lower
 *                    and upper indices of the ghost zone-centers/nodes
 *                    or edges at which data values should be assigned.
 * \param [in] side   specifies on which side boundary conditions need 
 *                    to be assigned. side can assume the following 
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
  double  *r, slp;

  if (side == X1_BEG){
    r = grid[IDIR].x;
    if (box->vpos == CENTER){
      BOX_LOOP(box,k,j,i){
        slp = r[i]/r[IBEG];
        d->Vc[RHO][k][j][i] = d->Vc[RHO][k][j][IBEG];
        d->Vc[VX1][k][j][i] = slp*d->Vc[VX1][k][j][IBEG];
        d->Vc[VX2][k][j][i] = slp*d->Vc[VX2][k][j][IBEG]; 
        d->Vc[PRS][k][j][i] = d->Vc[PRS][k][j][IBEG];
        d->Vc[BX1][k][j][i] = d->Vc[BX1][k][j][IBEG];
        d->Vc[BX2][k][j][i] = d->Vc[BX2][k][j][IBEG];
      }
    }else if (box->vpos == X2FACE){
      #ifdef STAGGERED_MHD
       BOX_LOOP(box,k,j,i) d->Vs[BX2s][k][j][i] = d->Vs[BX2s][k][j][IBEG];
      #endif
    }
  }
}

