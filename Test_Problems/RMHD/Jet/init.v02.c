/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Relativistic magnetized jet in axisymmetric coordinates.

  \b Description.\n
  The jet problem is set up in axisymmetric cylindrical coordinates
  \f$ (R,z) \f$.
  We label ambient and jet values with the suffix "a" and "j", respectively.
  The ambient medium is at rest and is characterized by uniform density
  and pressure, \f$ \rho_a \f$ and \f$ p_a \f$.
  The beam enters from the lower-z boundary from a circular nozzle of
  radius \f$ R_j \f$ and carries a constant poloidal field \f$ B_z \f$
  and a (radially-varying) toroidal component \f$ B_\phi(R) \f$.
  Flow variables are prescribed as
  \f[ \left\{\begin{array}{lcl}
     \rho(R)   &=& \rho_j  \\ \noalign{\medskip}
     v_z(R)    &=& v_j     \\ \noalign{\medskip}
     B_\phi(R) &=& \DS \left\{\begin{array}{ll}
             -\gamma B_m R/a & \quad{\rm for}\quad R<a \\ \noalign{\medskip}
             -\gamma B_m a/R & \quad{\rm for}\quad R>a
             \end{array}\right. \\ \noalign{\medskip}
     B_z(R)    &=& B_{z0}\quad\mathrm{(const)}
  \end{array}\right. \f]
  Here \c a is the magnetization radius while \f$ v_R = B_R = 0 \f$.
  These profiles are similar to the ones used by Tesileanu et al. (2008).
  The pressure distribution is found by solving the radial momentum
  balance between thermal, centrifugal and magnetic forces:
  \f[
    \frac{dp}{dR} = \frac{\rho v_\phi^2}{R} -
                    \frac{1}{2}\left[\frac{1}{R^2}\frac{d(R^2B^2_\phi)}{dR}
                                    +\frac{dB_z^2}{dR}\right]
  \f]
  Neglecting rotation and assuming \c Bz to be constant the solution of
  radial momentum balance becomes
  \f[
    p(R) = p_a + B_m^2\left[1-\min\left(\frac{R^2}{a^2},1\right)\right]
  \f]
  and the jet on-axis pressure increases for increasing toroidal field:
  \f[
    p(R=0) \equiv p_j = p_a + B_m^2
  \f]
  where \f$p_a\f$ is the ambient pressure.

  \b Normalization.\n
  Since the MHD equations are scale-invariant we have the freedom to
  specify a reference length, density and velocity. Here we choose
  - Length: jet radius \f$R_j=1\f$;
  - Density: ambient density \f$\rho_a=1\f$;
  - Velocity:
    - adiabtic setup (no cooling):
      ambient sound speed: \f$c_a = \Gamma p_a/\rho_a = 1\f$
      (from which it follows \f$p_a = 1/\Gamma\f$).
    - radiative setup (with cooling): 1 Km/s (set by \c UNIT_VELOCITY).
      In this case, the ambient pressure is computed from the ambient
      temperature <tt> Ta = 2500 K </tt>.

  In this way the number of parameters is reduced to 4:
  -# <tt>g_inputParam[ETA]</tt>: density contrast
     \f$ \eta = \rho_j/\rho_a
         \qquad\Longrightarrow\qquad \rho_j = \eta
     \f$
  -# <tt>g_inputParam[JET_VEL]</tt>:  Jet velocity  \f$ v_j \f$.
     This is also the Mach number of the beam with respect to the ambient
     in the adiabatic setup;
  -# <tt>g_inpurParam[SIGMA_Z]</tt>: poloidal magnetization;
  -# <tt>g_inpurParam[SIGMA_PHI]</tt>: toroidal magnetization;
  Magnetization are defined as follows:
  \f[ \left\{\begin{array}{ll}
       \sigma_z    &= \DS \frac{B_z^2}{2p_a}      \\ \noalign{\medskip}
       \sigma_\phi &= \DS \frac{<B_\phi^2>}{2p_a}
     \end{array}\right.
   \qquad\Longrightarrow\qquad
     \left\{\begin{array}{ll}
        B_z^2 &= \DS 2\sigma_zp_a       \\ \noalign{\medskip}
        B_m^2 &= \DS \frac{2\sigma_\phi}{a^2(1/2-2\log a)}p_a
      \end{array}\right.
  \f]
  Here the average value of \f$B^2_\phi\f$ is simply found as
  \f[
     <B^2_\phi> = \frac{\int_0^1 B^2_\phi R\,dR}{\int_0^1 R\,dR}
                = B_m^2a^2\left(\frac{1}{2} - 2\log a\right)
  \f]
  The following MAPLE code can be used to verify the solution:
  \code
   restart;   # Solve radial equilibrium balance equation
   assume (a < 1, a > 0);
   B[phi] := -B[m]*R/a;
   ode    := diff(p(R),R) = -1/(2*R^2)*diff(R^2*B[phi]^2,R);
   dsolve({ode,p(a)=pa});   # only for R < a
  \endcode

  When cooling is enabled, two additional parameters controlling the
  amplitude and frequency of perturbation can be used:
  -# <tt>g_inputParam[PERT_AMPLITUDE]</tt>: perturbation amplitude (in
     units of jet velocity);
  -# <tt>g_inputParam[PERT_PERIOD]</tt>: perturbation period (in years).

  - configurations #01-06 are adiabatic (no cooling);
  - configuration #07, #08 and #09 employs \c SNEq, \c MINEq and \c H2_COOL
    radiative cooling;

  The following image show the (log) density map at the end of
  simulation for setup #01.

  \image html mhd_jet.01.jpg "Density map at the end of the computation for configuration #01"

  \author A. Mignone (mignone@ph.unito.it)
  \date   Sept 12, 2014

  \b References
     - 
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

static void   GetJetValues (double x1, double x2, double x3, double *vj);
static double pmin, a;

/* ********************************************************************* */
void Init (double *v, double x1, double x2, double x3)
/*
 *
 *
 *
 *********************************************************************** */
{
  int    nv;
  double r, z, p0, vjet[256], vamb[256];

  g_gamma = 5./3.;
  GetJetValues (x1, x2, x3, vjet);

/* -- set ambient values -- */

  v[RHO] = g_inputParam[RHOA];
  EXPAND(v[VX1] = 0.0;   ,
         v[VX2] = 0.0;   ,
         v[VX3] = 0.0;)
  v[PRS] = pmin;
  #if PHYSICS == RMHD
   EXPAND(v[BX1] = 0.0;       ,
          v[BX2] = vjet[BX2];  ,
          v[BX3] = 0.0;)
   v[AX1] = 0.0;
   v[AX2] = 0.0;
   v[AX3] = vjet[AX3];
  #endif
  #if NTRACER > 0
   v[TRC] = 0.0;
  #endif
  #ifdef PSI_GLM 
   v[PSI_GLM] = 0.0;
  #endif

  g_smallPressure = 1.e-5;
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
/* 
 *
 *
 *********************************************************************** */
{
  int   i, j, k, nv;
  double *x1, *x2, *x3;
  double ***bxs, ***bys, ***bzs;
  double r, phi, prof, vjet[256], vout[NVAR];

  #ifdef STAGGERED_MHD
   D_EXPAND(bxs = d->Vs[BX1s];  , 
            bys = d->Vs[BX2s];  , 
            bzs = d->Vs[BX3s];)
  #endif

  x1  = grid[IDIR].xgc;  
  x2  = grid[JDIR].xgc;  
  x3  = grid[KDIR].xgc;

  if (side == X2_BEG){  /* -- X2_BEG boundary -- */
    if (box->vpos == CENTER){    /* -- cell-centered boundary conditions -- */
      BOX_LOOP(box,k,j,i){
 
        GetJetValues (x1[i], x2[j], x3[k], vjet); /* Jet Values */
        VAR_LOOP(nv) vout[nv] = d->Vc[nv][k][2*JBEG - j - 1][i]; /* Ambient */

        vout[VX2] *= -1.0; 
        EXPAND(vout[BX1] *= -1.0;  ,  
                                ;  , 
               vout[BX3] *= -1.0;)
        #ifdef PSI_GLM
         vjet[PSI_GLM]  = d->Vc[PSI_GLM][k][JBEG][i] - d->Vc[BX2][k][JBEG][i];
         vout[PSI_GLM] *= -1.0;
        #endif

        prof = (fabs(x1[i]) <= 1.0); // Profile(x1[i],x2[j],x3[k]);
        VAR_LOOP(nv) d->Vc[nv][k][j][i] = vout[nv] - (vout[nv] - vjet[nv])*prof;
      }

    }else if (box->vpos == X1FACE){  /* -- staggered fields -- */
      #ifdef STAGGERED_MHD
       x1 = grid[IDIR].xr;
       BOX_LOOP(box,k,j,i){
          vout[BX1] = -bxs[k][2*JBEG - j - 1][i];
          prof = (fabs(x1[i]) <= 1.0); // Profile(x1[i],x2[j],x3[k]);
          bxs[k][j][i] = vout[BX1] - (vout[BX1] - vjet[BX1])*prof;
       }
      #endif
    }
  }
}

/* **************************************************************** */
void GetJetValues (double x1, double x2, double x3, double *vj)
/*
 *
 *
 *
 * 
 ****************************************************************** */
{
  static int  first_call = 1;
  double bm, b2_av, bphi, lor;
  double r, x, scrh;

  a   = 0.5;
  r   = x1;
  lor = g_inputParam[LORENTZ];

  if (fabs(r) < 1.e-9) r = 1.e-9;

  x = r/a;
  vj[RHO] = g_inputParam[RHOJ];

  EXPAND(vj[VX1] = 0.0;                        ,
         vj[VX2] = sqrt(1.0 - 1.0/(lor*lor));  ,
         vj[VX3] = 0.0;)

  scrh = g_inputParam[MACH]/vj[VX2];
  pmin = vj[RHO]/(g_gamma*scrh*scrh - g_gamma/(g_gamma - 1.0));

  bm  = 4.0*pmin*g_inputParam[SIGMA_TOR];
  bm /= a*a*(1.0 - 4.0*log(a));
  bm  = sqrt(bm);

  scrh    = MIN(x*x, 1.0);
  vj[PRS] = pmin + bm*bm*(1.0 - scrh);

  bphi  = bm*(fabs(x) < 1.0 ? x: 1.0/x);
  bphi *= (r <= 1.0 ? 1.0:0.0);

  #if GEOMETRY == CYLINDRICAL
   EXPAND(vj[iBR]   = 0.0;                                     ,
          vj[iBZ]   = sqrt(2.0*g_inputParam[SIGMA_POL]*pmin);  ,
          vj[iBPHI] = lor*bphi;)
   vj[AX1] = 0.0;
   vj[AX2] = 0.0;
   vj[AX3] = 0.5*r*vj[iBZ];
  #endif

  if (first_call){
    print1 ("pmin  = %f\n",pmin);
    first_call = 0;
  }

/* -- set tracer value -- */

  #if NTRACER > 0
   vj[TRC] = 1.0;
  #endif
  #ifdef PSI_GLM 
   vj[PSI_GLM] = 0.0;
  #endif
}
