/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Contains some implementation of the PatchPluto class.

   PatchTools.cpp contains a number of function used by PatchUnsplit.
   Most of them are simple wrappers to other functions employed by the 
   static version of the code.

  \authors C. Zanni   (zanni@oato.inaf.it)\n
           A. Mignone (mignone@ph.unito.it)
  \date   Sep 20, 2012
*/
/* ///////////////////////////////////////////////////////////////////// */

#include <cstdio>
#include <string>
using std::string;

#include "PatchPluto.H"
#include "LoHiSide.H"

/* ********************************************************************* */
void PatchPluto::saveFluxes (const State_1D *state, int beg, int end, 
                             Grid *grid)
/*! 
 *  Rebuild fluxes in a way suitable for AMR operation
 *  by adding pressure and multiplying by area.
 *
 *********************************************************************** */
{
  int  i, nv;
  double **f, *p, r, area;

  f = state->flux;
  p = state->press;

  for (i = beg; i <= end; i++) f[i][MXn] += p[i];

  #if (GEOMETRY == CARTESIAN) && (CH_SPACEDIM > 1)
    if ((g_dir == IDIR) && (g_stretch_fact != 1.)) {
      for (i = beg; i <= end; i++) {
        VAR_LOOP(nv) f[i][nv] *= g_stretch_fact;
      }
    }
   #if (CH_SPACEDIM == 3)
    if ((g_dir == JDIR) && (g_x3stretch != 1.)) {
      for (i = beg; i <= end; i++) {
        VAR_LOOP(nv) f[i][nv] *= g_x3stretch;
      }
    }
    if ((g_dir == KDIR) && (g_x2stretch != 1.)) {
      for (i = beg; i <= end; i++) {
        VAR_LOOP(nv) f[i][nv] *= g_x2stretch;
      }
    }   
   #endif
  #endif

  #if GEOMETRY == CYLINDRICAL
   if (g_dir == IDIR){
     for (i = beg; i <= end; i++) {
     for (nv = 0; nv < NVAR; nv++) {
       f[i][nv] *= grid[IDIR].A[i];
       #if CH_SPACEDIM > 1
        f[i][nv] *= g_x2stretch;
       #endif
     }}
   }else{
     area = fabs(grid[IDIR].x[g_i]);
     for (i = beg; i <= end; i++) {
     for (nv = 0; nv < NVAR; nv++) {
       f[i][nv] *= area;
     }}
   }
  #endif

  #if GEOMETRY == SPHERICAL
   if (g_dir == IDIR){
    #if CH_SPACEDIM > 1
     area = grid[JDIR].dV[g_j]/m_dx;
     #if CH_SPACEDIM == 3
      area *= g_x3stretch;
     #endif 
    #endif
     for (i = beg; i <= end; i++) {
      for (nv = 0; nv < NVAR; nv++) {
        f[i][nv] *= grid[IDIR].A[i];
       #if CH_SPACEDIM > 1
        f[i][nv] *= area;
       #endif 
      }
      #if (COMPONENTS == 3) && (CHOMBO_EN_SWITCH == YES)
       f[i][iMPHI] *= grid[IDIR].xr[i]*sin(grid[JDIR].x[g_j]);
      #endif
     }
   }
   if (g_dir == JDIR){
     area = fabs(grid[IDIR].x[g_i]);
     #if CHOMBO_LOGR == YES
      area *= grid[IDIR].dx[g_i]/m_dx;
     #endif
     #if CH_SPACEDIM == 3
      area *= g_x3stretch;
     #endif
     for (i = beg; i <= end; i++) {
      for (nv = 0; nv < NVAR; nv++) {
        f[i][nv] *= grid[JDIR].A[i]*area;
      }
      #if (COMPONENTS == 3) && (CHOMBO_EN_SWITCH == YES)
       f[i][iMPHI] *= grid[IDIR].x[g_i]*sin(grid[JDIR].xr[i]);
      #endif
     }
   }
   if (g_dir == KDIR){
     area = g_x2stretch*fabs(grid[IDIR].x[g_i]);
    #if CHOMBO_LOGR == YES 
     area *= grid[IDIR].dx[g_i]/m_dx;
    #endif
     for (i = beg; i <= end; i++) {
      for (nv = 0; nv < NVAR; nv++) {
        f[i][nv] *= area;
      }
      #if (COMPONENTS == 3) && (CHOMBO_EN_SWITCH == YES)
       f[i][iMPHI] *= grid[IDIR].x[g_i]*sin(grid[JDIR].x[g_j]);
      #endif
     }
   }
  #endif

  #if GEOMETRY == POLAR
   if (g_dir == IDIR){
     for (i = beg; i <= end; i++) {
      for (nv = 0; nv < NVAR; nv++) {
        f[i][nv] *= grid[IDIR].A[i];
       #if CH_SPACEDIM > 1
        f[i][nv] *= g_x2stretch;
       #endif
       #if CH_SPACEDIM == 3
        f[i][nv] *= g_x3stretch;
       #endif
      }
      #if (COMPONENTS > 1) && (CHOMBO_EN_SWITCH == YES)
       f[i][iMPHI] *= grid[IDIR].xr[i];
      #endif
     }
   }
   if (g_dir == JDIR) {
     area = g_x3stretch;
     #if CHOMBO_LOGR == YES
      area *= grid[IDIR].dx[g_i]/m_dx;
     #endif
     if (area != 1.) {
      for (i = beg; i <= end; i++) {
      for (nv = 0; nv < NVAR; nv++) {
        f[i][nv] *= area;
     }}}
     #if (COMPONENTS > 1) && (CHOMBO_EN_SWITCH == YES)
      for (i = beg; i <= end; i++) f[i][iMPHI] *= grid[IDIR].x[g_i];
     #endif
   }
   if (g_dir == KDIR) {
     area = g_x2stretch*fabs(grid[IDIR].x[g_i]);
    #if CHOMBO_LOGR == YES
     area *= grid[IDIR].dx[g_i]/m_dx;
    #endif
     for (i = beg; i <= end; i++) {
      for (nv = 0; nv < NVAR; nv++) {
        f[i][nv] *= area;
      }
      #if (COMPONENTS > 1) && (CHOMBO_EN_SWITCH == YES)
       f[i][iMPHI] *= grid[IDIR].x[g_i];
      #endif
     }
   }
  #endif
}
#if DIMENSIONS > 1
/* ********************************************************************* */
void PatchPluto::getPrimitiveVars (Data_Arr U, Data *d, Grid *grid)
/*!
 * - Recover primitive variables from the input conservative array \c U
 * - Set physical boundary conditions and convert
 *
 *********************************************************************** */
{
  int i, j, k, nv, dir;
  int nx, ny, nz;
  int ibeg, jbeg, kbeg;
  int iend, jend, kend;
  int lft_side[3], rgt_side[3];
  int err;
  double cylr;

  nx = grid[IDIR].np_tot;
  ny = grid[JDIR].np_tot;
  nz = grid[KDIR].np_tot;

/* ------------------------------------------------------- 
     check whether the patch touches a physical boundary
   ------------------------------------------------------- */

  for (dir = 0; dir < DIMENSIONS; dir++){
    lft_side[dir] = (grid[dir - IDIR].lbound != 0);
    rgt_side[dir] = (grid[dir - IDIR].rbound != 0);
  }

/* ---------------------------------------------------
     Extract the portion of the domain where U 
     is defined (i.e. NOT in the physical boundary).
   --------------------------------------------------- */

  ibeg = 0; iend = nx - 1;
  jbeg = 0; jend = ny - 1;
  kbeg = 0; kend = nz - 1;

/* -------------------------------------------------
     exclude physical boundaries since the 
     conservative vector U is not yet defined.    
   ------------------------------------------------- */

  D_EXPAND(if (lft_side[IDIR]) ibeg = IBEG;  ,
           if (lft_side[JDIR]) jbeg = JBEG;  ,
           if (lft_side[KDIR]) kbeg = KBEG;)

  D_EXPAND(if (rgt_side[IDIR]) iend = IEND;  ,
           if (rgt_side[JDIR]) jend = JEND;  ,
           if (rgt_side[KDIR]) kend = KEND;)

/* ---------------------------------------------------
            recover total energy 
   --------------------------------------------------- */

  #if CHOMBO_EN_SWITCH == YES
   totEnergySwitch (U, ibeg, iend, jbeg, jend, kbeg, kend, +1);
  #endif

/* ----------------------------------------------
     convert from conservative to primitive 
   ---------------------------------------------- */

  convertConsToPrim(U, d->Vc, ibeg, jbeg, kbeg, iend, jend, kend, grid);

/* ----------------------------------------------
         Assign boundary conditions 
   ---------------------------------------------- */
    
  Boundary (d, ALL_DIR, grid);

/* --------------------------------------------------------------
    set boundary values on conservative variables by converting 
    primitive variables in the corresponding regions.
   -------------------------------------------------------------- */

  #if INTERNAL_BOUNDARY == YES
   convertPrimToCons(d->Vc, U, 0, 0, 0, NX1_TOT-1, NX2_TOT-1, NX3_TOT-1, grid);
  #else
 
   if (lft_side[IDIR]) convertPrimToCons(d->Vc, U, 0, 0, 0, 
                                         IBEG-1, NX2_TOT-1, NX3_TOT-1, grid);
  
   if (lft_side[JDIR]) convertPrimToCons(d->Vc, U, 0, 0, 0,
                                         NX1_TOT-1, JBEG-1, NX3_TOT-1, grid);
    
   if (lft_side[KDIR]) convertPrimToCons(d->Vc, U, 0, 0, 0,
                                         NX1_TOT-1, NX2_TOT-1, KBEG-1, grid);

   if (rgt_side[IDIR]) convertPrimToCons(d->Vc, U, IEND+1, 0, 0,
                                         NX1_TOT-1, NX2_TOT-1, NX3_TOT-1, grid);
  
   if (rgt_side[JDIR]) convertPrimToCons(d->Vc, U, 0, JEND+1, 0,
                                         NX1_TOT-1, NX2_TOT-1, NX3_TOT-1, grid);
    
   if (rgt_side[KDIR]) convertPrimToCons(d->Vc, U, 0, 0, KEND+1,
                                         NX1_TOT-1, NX2_TOT-1, NX3_TOT-1, grid);
  #endif
}
#endif /* DIMENSIONS > 1 */

/* ************************************************************ */
void PatchPluto::showPatch (Grid *grid)
/*
 *
 *
 *
 ************************************************************** */
{
  char pb[4]="*";  /* -- physical boundary -- */
  char ib[4]="o";  /* -- internal boundary -- */

  Grid *Gx, *Gy, *Gz;

  D_EXPAND(Gx = grid + IDIR;  ,
           Gy = grid + JDIR;  ,
           Gz = grid + KDIR;)

  print ("+-----------------------------------------------------------+\n");
  print ("| Level = %d \n", grid[IDIR].level);
  print ("| ib,ie = [%d, %d], xb,xe = %s[%f, %f]%s\n", 
             IBEG, IEND, (Gx->lbound != 0 ? pb:ib), 
              Gx->xr[IBEG-1], Gx->xr[IEND],(Gx->rbound != 0 ? pb:ib));
  print ("| jb,je = [%d, %d], yb,ye = %s[%f, %f]%s\n", 
            JBEG, JEND, (Gy->lbound != 0 ? pb:ib), 
            Gy->xr[JBEG-1], Gy->xr[JEND], (Gy->rbound != 0 ? pb:ib));
  print ("+-----------------------------------------------------------+\n");

}

/* ************************************************************** */
void PatchPluto::convertConsToPrim (Data_Arr U, Data_Arr V, 
                               int ibeg, int jbeg, int kbeg,
                               int iend, int jend, int kend, Grid *grid) 
/*!
 *  Implements as simple wrapper to 1D functions in order to 
 *  convert U (conservative) into d->Vc (primitive).
 *
 *  Furthermore: during the first step (g_intStage = 1, this 
 *  function is being called from getPrimitiveVars only), 
 *  re-normalize ions when using the MINEq chemical network.
 * 
 **************************************************************** */
{
  int i, j, k, nv, err;
  static double **u, **v;
  static unsigned char *flag;
  Index indx;

  if (u == NULL){
    u = ARRAY_2D(NMAX_POINT, NVAR, double);
    v = ARRAY_2D(NMAX_POINT, NVAR, double);
    flag = ARRAY_1D(NMAX_POINT, unsigned char);
  }

  g_dir = IDIR;
  SetIndexes(&indx, grid);
  for (k = kbeg; k <= kend; k++){ g_k = k;
  for (j = jbeg; j <= jend; j++){ g_j = j;
    for (i = ibeg; i <= iend; i++){
    for (nv = NVAR; nv--;   ){
      u[i][nv] = U[nv][k][j][i];
    }}

    #if COOLING == MINEq || COOLING == H2_COOL
     if (g_intStage == 1) {
       for (i = ibeg; i <= iend; i++) NormalizeIons(u[i]);
     }
    #endif
 
    err = ConsToPrim (u, v, ibeg, iend, flag);

    if (err != 0) {
      #if WARNING_MESSAGES == YES
       print (">> from convertConsToPrim, t = %f, ", g_time);
       for (i = ibeg; i <= iend; i++){ 
       if (flag[i] != 0) {
         #if DIMENSIONS == 2
          print ("zone [%d %d], ",i,j);
         #elif DIMENSIONS == 3
          print ("zone [%d %d ], ",i,j,k);
         #endif
       }}
       print("\n");
       showPatch(grid);
      #endif
    }

    for (i = ibeg; i <= iend; i++){
    for (nv = NVAR; nv--;   ){
      V[nv][k][j][i] = v[i][nv];
    }}
  }}
}

/* ************************************************************** */
void PatchPluto::convertPrimToCons (Data_Arr V, Data_Arr U,
                               int ibeg, int jbeg, int kbeg,
                               int iend, int jend, int kend, Grid *grid) 
/*
 *
 *  Convert d->Vc (primitive) into U (conservative).
 *  This is a simple wrapper to 1D functions.
 *
 *
 **************************************************************** */
{
  int i, j, k, nv, err;
  static double **u, **v;
  static unsigned char *flag;
  Index indx;

  if (u == NULL){
    u = ARRAY_2D(NMAX_POINT, NVAR, double);
    v = ARRAY_2D(NMAX_POINT, NVAR, double);
    flag = ARRAY_1D(NMAX_POINT, unsigned char);
  }

  g_dir = IDIR;
  SetIndexes(&indx, grid);
  for (k = kbeg; k <= kend; k++){ g_k = k;
  for (j = jbeg; j <= jend; j++){ g_j = j;
    for (i = ibeg; i <= iend; i++){
    for (nv = NVAR; nv--;   ){
      v[i][nv] = V[nv][k][j][i];
    }}
    PrimToCons(v, u, ibeg, iend);
    for (i = ibeg; i <= iend; i++){
    for (nv = NVAR; nv--;   ){
      U[nv][k][j][i] = u[i][nv];
    }}
  }}
}
/* ********************************************************************* */
void PatchPluto::convertFArrayBox(FArrayBox&  U)
/*!
 *  Convert a conservative array to primitive for 
 *  plotfile data. 
 *  Called from AMRLevelPluto::writePlotLevel
 *
 *
 *********************************************************************** */
{
  int ibeg, jbeg, kbeg;
  int iend, jend, kend;
  int i,j,k, nv;
  double ***UU[NVAR];
  static unsigned char *flag;
  static double **u, **v;

  if (u == NULL){
    u    = ARRAY_2D(NMAX_POINT, NVAR, double);
    v    = ARRAY_2D(NMAX_POINT, NVAR, double);
    flag = ARRAY_1D(NMAX_POINT, unsigned char);
  }

  jbeg = jend = kbeg = kend = 0;
  D_EXPAND(ibeg = U.loVect()[IDIR]; 
           iend = U.hiVect()[IDIR]; ,
           jbeg = U.loVect()[JDIR];   
           jend = U.hiVect()[JDIR]; ,
           kbeg = U.loVect()[KDIR];   
           kend = U.hiVect()[KDIR]; );

  for (nv=0; nv<NVAR; nv++) {
    UU[nv] = ArrayBoxMap(kbeg,kend,jbeg,jend,ibeg,iend,U.dataPtr(nv));
  }

/* ----------------------------------------------------------------
    HOTIFX: set g_intStage to a negative value so that ConsToPrim
    does not have to call CheckZone when writing HDF5 data 
    to disk since d->flag is not defined or does not
    overlap with U (which is an FArrayBox)
   ---------------------------------------------------------------- */

  g_intStage = -1;

/* --------------------------------------------------------
        Recover total energy if necessary.
   -------------------------------------------------------- */

  #if CHOMBO_EN_SWITCH == YES
   totEnergySwitch (UU, ibeg+IOFFSET, iend-IOFFSET, 
                        jbeg+JOFFSET, jend-JOFFSET, 
                        kbeg+KOFFSET, kend-KOFFSET, +1);
  #endif

/* --------------------------------------------------------
    Conversion is done in the interior points only. 
    For this reason we exclude from (ibeg, iend, ...) one 
    boundary zone in each direction.
   -------------------------------------------------------- */

  for (k = kbeg + KOFFSET; k <= kend - KOFFSET; k++){
  for (j = jbeg + JOFFSET; j <= jend - JOFFSET; j++){
    for (i = ibeg + IOFFSET; i <= iend - IOFFSET; i++){
    for (nv=0; nv < NVAR; nv++){
      u[i-ibeg][nv] = UU[nv][k][j][i];
    }}
    ConsToPrim (u, v, IOFFSET, iend-ibeg-IOFFSET, flag);
    for (i = ibeg; i <= iend; i++){
    for (nv=0 ; nv < NVAR; nv++){
      UU[nv][k][j][i] = v[i-ibeg][nv];
    }}
  }}

  for (nv = 0; nv < NVAR; nv++) 
    FreeArrayBoxMap(UU[nv], kbeg, kend, jbeg, jend, ibeg, iend);
}

/* ********************************************************************* */
void PatchPluto::totEnergySwitch(Data_Arr UU, int ibeg, int iend, 
                 int jbeg, int jend, int kbeg, int kend, int s)
/*! 
 *  Perform a number of opeations that involve replacing energy 
*   with entropy or viceversa.
 *  In particular, 
 *
 *  - if s =  1, convert {D, m, B, sigma_c} -> {D, m, B, E}
 *  - if s = -1, convert {D, m, B,       E} -> {D, m, B, sigma_c}
 *  - if s =  0, there're two possible entropy values: one is obtained from 
 *               total energy sigma_c(E) (stored in UU[ENG]) and the other one 
 *               is the actual entropy variable sigma_c updated by enabling 
 *               ENTROPY_SWITCH.
 *               if FLAG_ENTROPY is true -> sigma_c(E) = sigma_c
 *               otherwise               -> sigma_c    = sigma_c(E)
 *
 *
 *********************************************************************** */
{
#if HAVE_ENERGY
  int i,j,k;
  double D, mx, my, mz, E, sigma_c;
  double Bx, By, Bz;
  double p, sigma, Kin;
  #if PHYSICS == RMHD || PHYSICS == RHD
   static double **v;
   Map_param par;

   if (v == NULL) v = ARRAY_2D(1, NVAR, double);
  #endif

 /* -------------------------------------------------------------
      s == 1: convert {D, m, B, sigma_c} -> {D, m, B, E}
    ------------------------------------------------------------- */

  if (s == 1){
    for (k = kbeg; k <= kend; k++){
    for (j = jbeg; j <= jend; j++){
    for (i = ibeg; i <= iend; i++){
      D       = UU[RHO][k][j][i];
      sigma_c = UU[ENG][k][j][i];
      EXPAND(mx = UU[MX1][k][j][i];  ,
             my = UU[MX2][k][j][i];  ,
             mz = UU[MX3][k][j][i];)

      #if PHYSICS == MHD || PHYSICS == RMHD
       EXPAND(Bx = UU[BX1][k][j][i];  ,
              By = UU[BX2][k][j][i];  ,
              Bz = UU[BX3][k][j][i];)
      #endif

      #if PHYSICS == RHD || PHYSICS == RMHD

       par.D       = D;
       par.sigma_c = sigma_c;
       par.m2 = EXPAND(mx*mx, + my*my, + mz*mz);
       #if PHYSICS == RMHD
        par.S  = EXPAND(mx*Bx, + my*By, + mz*Bz);
        par.B2 = EXPAND(Bx*Bx, + By*By, + Bz*Bz); 
        par.S2 = par.S*par.S;
       #endif
       EntropySolve(&par);
       UU[ENG][k][j][i] = par.E;

      #else

       Kin  = EXPAND(mx*mx, + my*my, + mz*mz);
       Kin /= D;
       #if PHYSICS == MHD
        Kin += EXPAND(Bx*Bx, + By*By, + Bz*Bz);
       #endif
       #if EOS == IDEAL
        p = sigma_c*pow(D,g_gamma-1.0);
        UU[ENG][k][j][i] = p/(g_gamma - 1.0) + 0.5*Kin;
       #else
        print1 ("! totEnergySwitch: CHOMBO_EN_SWITCH only works with IDEAL EOS for MHD/HD\n");
        QUIT_PLUTO(1);
       #endif
      #endif
    }}}
  }
 
 /* -------------------------------------------------------------
      s == -1: convert {D, m, B, E} -> {D, m, B, sigma_c}
    ------------------------------------------------------------- */

  if (s == -1){
    for (k = kbeg; k <= kend; k++){
    for (j = jbeg; j <= jend; j++){
    for (i = ibeg; i <= iend; i++){
      D = UU[RHO][k][j][i];
      E = UU[ENG][k][j][i];
      EXPAND(mx = UU[MX1][k][j][i];  ,
             my = UU[MX2][k][j][i];  ,
             mz = UU[MX3][k][j][i];)

      #if PHYSICS == MHD || PHYSICS == RMHD
       EXPAND(Bx = UU[BX1][k][j][i];  ,
              By = UU[BX2][k][j][i];  ,
              Bz = UU[BX3][k][j][i];)
      #endif

      #if PHYSICS == RHD || PHYSICS == RMHD

       par.D  = D;
       par.E  = E;
       par.m2 = EXPAND(mx*mx, + my*my, + mz*mz);
       #if PHYSICS == RMHD
        par.S  = EXPAND(mx*Bx, + my*By, + mz*Bz);
        par.B2 = EXPAND(Bx*Bx, + By*By, + Bz*Bz); 
        par.S2 = par.S*par.S;
       #endif
       EnergySolve(&par);
       v[0][RHO] = par.rho;  /* par.D/par.lor */
       v[0][PRS] = par.prs;
       Entropy (v, &sigma, 0, 0);
       UU[ENG][k][j][i] = par.D*sigma;

      #else

       Kin  = EXPAND(mx*mx, + my*my, + mz*mz);
       Kin /= D;
       #if PHYSICS == MHD
        Kin += EXPAND(Bx*Bx, + By*By, + Bz*Bz);
       #endif
       #if EOS == IDEAL
        p = (g_gamma - 1.0)*(UU[ENG][k][j][i] - 0.5*Kin);
        UU[ENG][k][j][i] = p/pow(D,g_gamma-1.0);
       #else
        print1 ("! totEnergySwitch: CHOMBO_EN_SWITCH only works with IDEAL EOS for MHD/HD\n");
        QUIT_PLUTO(1);
       #endif

      #endif
    }}}
  }

 /* ----------------------------------------------------------------
      s == 0: replace sigma_c(E) with sigma_c if FLAG_ENTROPY is on
              otherwise sigma_c = sigma_c(E)
              IMPORTANT: call this function with s=0 only after
              energy has been converted into entropy!
    ---------------------------------------------------------------- */

  #if ENTROPY_SWITCH == YES
   #if CHOMBO_EN_SWITCH == NO
    #error CHOMBO_EN_SWITCH should be turned on when using ENTROPY_SWITCH
   #endif
   
   if (s == 0){
    #if (COOLING == NO) 
     g_dir = IDIR; 
     for (k = kbeg; k <= kend; k++){ g_k = k;
     for (j = jbeg; j <= jend; j++){ g_j = j;
     for (i = ibeg; i <= iend; i++){ g_i = i;
       if (CheckZone(i, FLAG_ENTROPY)){
         UU[ENG][k][j][i] = UU[ENTR][k][j][i];
       }else{
         UU[ENTR][k][j][i] = UU[ENG][k][j][i];
       }      
     }}}
    #else
     for (k = kbeg; k <= kend; k++){
     for (j = jbeg; j <= jend; j++){
     for (i = ibeg; i <= iend; i++){
       UU[ENTR][k][j][i] = UU[ENG][k][j][i];
     }}}
    #endif
   }
  #endif
#endif
}
