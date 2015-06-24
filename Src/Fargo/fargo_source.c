#include "pluto.h"

/* ********************************************************************* */
void FARGO_ADD_SOURCE(const State_1D *state, int beg, int end, 
                      double dt, Grid *grid)
/* 
 * PURPOSE
 *
 *   Add source terms to the azimuthal component of the momentum
 *   equation.
 * 
 *********************************************************************** */
{
  int    i;
  double w, dw, vphi;
  double **rhs, *u, *v, **wA;
  double *r, *ct, *dx;

  wA  = FARGO_GetVelocity();
  rhs = state->rhs;
  dx  = grid[g_dir].dx;

  #if GEOMETRY == CARTESIAN

   if (g_dir == IDIR){

     for (i = beg; i <= end; i++){
       v  = state->vh[i]; 
       dw = 0.5*(wA[g_k][i+1] - wA[g_k][i-1])/dx[i];

       rhs[i][MX2] -= dt*u[MX1]*dw;
       #if HAVE_ENERGY
        rhs[i][ENG] -= dt*v[VX2]*u[MX1]*dw;
        #if PHYSICS == MHD
         rhs[i][ENG] += dt*v[BX2]*v[BX1]*dw;
        #endif
       #endif
     }

   } else if (g_dir == KDIR){

     for (i = beg; i <= end; i++){
       v  = state->vh[i]; u = state->uh[i];
       dw = 0.5*(wA[i+1][g_i] - wA[i-1][g_i])/dx[i];
       
       rhs[i][MX2] -= dt*u[MX3]*dw;
       #if HAVE_ENERGY
        rhs[i][ENG] -= dt*v[VX2]*u[MX3]*dw;
        #if PHYSICS == MHD
         rhs[i][ENG] += dt*v[BX2]*v[BX3]*dw;
        #endif
       #endif
     }
   } 
   
  #elif GEOMETRY == POLAR

   r = grid[IDIR].x;
   if (g_dir == IDIR){

     for (i = beg; i <= end; i++){
       v  = state->vh[i]; u = state->uh[i];
       w  = wA[g_k][i];
       dw = 0.5*(wA[g_k][i+1] - wA[g_k][i-1])/dx[i];

       rhs[i][iMPHI] -= dt*u[iMR]*(dw + w/r[i]);
       #if HAVE_ENERGY
        vphi = v[iVPHI] + w;
        rhs[i][ENG] += dt*u[iMR]*(vphi*w/r[i] - v[iVPHI]*dw);
        #if PHYSICS == MHD
         rhs[i][ENG] -= dt*v[iBR]*v[iBPHI]*(w/r[i] - dw);
        #endif
       #endif
     }

   }else if (g_dir == KDIR){

     for (i = beg; i <= end; i++){
       v  = state->vh[i]; u = state->uh[i];
       w  = wA[i][g_i];
       dw = 0.5*(wA[i+1][g_i] - wA[i-1][g_i])/dx[i];
       
       rhs[i][iMPHI] -= dt*u[MX3]*dw;
       #if HAVE_ENERGY
        rhs[i][ENG] -= dt*v[VX2]*u[MX3]*dw;
        #if PHYSICS == MHD
         rhs[i][ENG] += dt*v[BX2]*v[BX3]*dw;
        #endif
       #endif
     }
   }

  #elif GEOMETRY == SPHERICAL

   r  = grid[IDIR].x;
   if (g_dir == IDIR){

     for (i = beg; i <= end; i++){
       v  = state->vh[i]; u = state->uh[i];
       w  = wA[g_j][i];
       dw = 0.5*(wA[g_j][i+1] - wA[g_j][i-1])/dx[i];

       rhs[i][iMPHI] -= dt*u[iMR]*(dw + w/r[i]);
       #if HAVE_ENERGY
        vphi = v[iVPHI] + w;
        rhs[i][ENG] += dt*u[iMR]*(vphi*w/r[i] - v[iVPHI]*dw);
        #if PHYSICS == MHD
         rhs[i][ENG] -= dt*v[iBR]*v[iBPHI]*(w/r[i] - dw); 
        #endif
       #endif
     }

   }else if (g_dir == JDIR){

     ct = grid[JDIR].ct;
     for (i = beg; i <= end; i++){
       v  = state->vh[i]; u = state->uh[i];
       w  = wA[i][g_i];
       dw = 0.5*(wA[i+1][g_i] - wA[i-1][g_i])/(r[g_i]*dx[i]);

       rhs[i][iMPHI] -= dt*u[iMTH]*(dw + ct[i]*w/r[g_i]);

       #if HAVE_ENERGY
        vphi = v[iVPHI] + w;
        rhs[i][ENG] += dt*u[iMTH]*(ct[i]*vphi*w/r[i] - v[iVPHI]*dw);
        #if PHYSICS == MHD
         rhs[i][ENG] -= dt*v[iBTH]*v[iBPHI]*(ct[i]*w/r[g_i] - dw); 
        #endif
       #endif
     }
   }
  #endif /* GEOMETRY == SPHERICAL */
}
