#include "pluto.h"

/* ********************************************************************* */
void Init (double *us, double x, double y, double z)
/*
 *
 *
 *
 *********************************************************************** */
{
  double scrh = 1.0;

  g_gamma = g_inputParam[GAMMA_EOS];
  if (x < 0.5){
    us[RHO] = g_inputParam[DN_L];
    us[VX1] = g_inputParam[VX_L];
    us[VX2] = g_inputParam[VY_L];
    us[PRS] = g_inputParam[PR_L];
  }else{
    us[RHO] = g_inputParam[DN_R];
    us[VX1] = g_inputParam[VX_R];
    us[VX2] = g_inputParam[VY_R];
    us[PRS] = g_inputParam[PR_R];
  }

  #if USE_FOUR_VELOCITY == YES
   scrh = EXPAND(us[VX1]*us[VX1], + us[VX2]*us[VX2], + 0.0);
   scrh = 1.0/sqrt(1.0 - scrh);
   us[VX1] *= scrh;
   us[VX2] *= scrh;
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
{ }

