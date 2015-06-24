/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Flag zones where pressure is updated using entropy.

  Flag cells where pressure has to be updated from the conserved
  entropy density and compute entropy at the beginning of each
  integration stage.
  
  The flagging strategy is based on two switches designed to detect 
  the presence of compressive motion or shock waves in the fluid:
  - \f$ \nabla\cdot\vec{v} < 0 \f$
  - \f$ |\nabla p|/p > 0.1 \f$
  
  By default, if at least one of the switch is \e FALSE, we flag the 
  computational zone by turning the ::FLAG_ENTROPY bit of the 
  \c d->flag array on.
  In this way, the default strategy is to evolve the entropy equation
  everywhere, except at shocks (where both switches are \e TRUE).

  This flag will be checked later in the ConsToPrim() functions in 
  order to recover pressure from the entropy density rather than 
  from the total energy density.
  The update process is:
 
  - start with a vector of primitive variables  <tt> {V,s} </tt>
    where \c s is the entropy;
  - set boundary condition on \c {V}; 
    compute \c {s} from \c {V};
    flag zones where entropy may be used (flag = 1);
  - evolve equations for one time step;
  - convert <tt> {U,S} </tt> to primitive:
   
        where (flag == 0) {
          p = p(E);
        }else{  
          p = p(S);
          correct E using new p;
        }
                         
  \b Reference
     - "Maintaining Pressure Positivity in Magnetohydrodynamics Simulations"
       Balsara \& Spicer, JCP (1999) 148, 133
 

  \author A. Mignone (mignone@ph.unito.it)
  \date   Aug 16, 2012
*/
/* ///////////////////////////////////////////////////////////////////// */
#include"pluto.h"

#if HAVE_ENERGY && (ENTROPY_SWITCH == YES)
/* ********************************************************************* */
void ComputeEntropy (const Data *d, Grid *grid)
/*!
 * Compute entropy as a primitive variable.
 *
 * \param [in,out]  d   pointer to Data structure
 * \param [in]    grid  pointer to an array of Grid structures.
 *
 * \return This function has no return value.
 *
 *********************************************************************** */
{
  int i, j, k;
  static double **v1d;

  if (v1d == NULL) v1d = ARRAY_2D(NMAX_POINT, NVAR, double);
  KTOT_LOOP(k) {
  JTOT_LOOP(j) {
    ITOT_LOOP(i) {
      v1d[i][RHO] = d->Vc[RHO][k][j][i];
      v1d[i][PRS] = d->Vc[PRS][k][j][i];
    }
    Entropy(v1d, d->Vc[ENTR][k][j], 0, NX1_TOT-1);
  }}
}

/* ********************************************************************* */
void EntropyOhmicHeating (const Data *d, Data_Arr UU, double dt, Grid *grid)
/*! 
 * Add Ohmic heating term to the conservative entropy equation when
 * RESISTIVE_MHD is set to YES.
 *
 * \return  This function has no return value.
 *********************************************************************** */
{
#if PHYSICS == MHD && (RESISTIVE_MHD == EXPLICIT) /* at the moment ... remove later ! */
  int    i,j,k,nv;
  double rho, rhog, gm1, vc[NVAR];
  double Jc[3], eta[3], J2eta;
  double ***Jx1, *x1, *dx1;
  double ***Jx2, *x2, *dx2;
  double ***Jx3, *x3, *dx3;
  Data_Arr V;

/* ---------------------------------------------
   1. Set pointer shortcuts 
   --------------------------------------------- */
  
  V   = d->Vc;
  Jx1 = d->J[IDIR]; x1 = grid[IDIR].x ; dx1 = grid[IDIR].dx;
  Jx2 = d->J[JDIR]; x2 = grid[JDIR].x ; dx2 = grid[JDIR].dx;
  Jx3 = d->J[KDIR]; x3 = grid[KDIR].x ; dx3 = grid[KDIR].dx;

  gm1 = g_gamma - 1.0;
  Jc[IDIR] = Jc[JDIR] = Jc[KDIR] = 0.0;
  
/* ---------------------------------------------
   2. Main spatial loop
   --------------------------------------------- */
  
  DOM_LOOP(k,j,i){

  /* ---------------------------------
     2a. compute currents at cell center
     --------------------------------- */

    #ifdef STAGGERED_MHD  /* Staggered MHD */

     #if COMPONENTS == 3
      Jc[IDIR] = AVERAGE_YZ(Jx1,k-1,j-1,i);
      Jc[JDIR] = AVERAGE_XZ(Jx2,k-1,j,i-1);
     #endif
     Jc[KDIR] = AVERAGE_XY(Jx3,k,j-1,i-1);

    #else               /* Cell-centered MHD */

     #if COMPONENTS == 3
      Jc[IDIR] = CDIFF_X2(V[BX3],k,j,i)/dx2[j] - CDIFF_X3(V[BX2],k,j,i)/dx3[k];
      Jc[JDIR] = CDIFF_X3(V[BX1],k,j,i)/dx3[k] - CDIFF_X1(V[BX3],k,j,i)/dx1[i];
     #endif
     Jc[KDIR] = CDIFF_X1(V[BX2],k,j,i)/dx1[i] - CDIFF_X2(V[BX1],k,j,i)/dx2[j];

     #if GEOMETRY != CARTESIAN
      print1 ("! EntropyOhmicHeating: only CT supported in this geometry.\n")
      QUIT_PLUTO(1);
     #endif

    #endif

  /* ----------------------------------------
     2b. compute resistivity at cell center.
     ---------------------------------------- */

    VAR_LOOP(nv) vc[nv] = V[nv][k][j][i];
    rho  = vc[RHO];
    rhog = pow(rho,-gm1);

    Resistive_eta (vc, x1[i], x2[j], x3[k], Jc, eta);
    J2eta = 0.0;
    for (nv = 0; nv < 3; nv++) J2eta += Jc[nv]*Jc[nv]*eta[nv];

  /* ----------------------------------------
     2c. Update conserved entropy
     ---------------------------------------- */
     
    UU[k][j][i][ENTR] += dt*rhog*gm1*J2eta;
  }
#endif
}
#endif
