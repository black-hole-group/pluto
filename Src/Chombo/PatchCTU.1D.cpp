#include <cstdio>
#include <string>
using std::string;

#include "PatchPluto.H"
#include "LoHiSide.H"

/* *********************************************************** */
void PatchPluto::updateSolution(FArrayBox&       a_U,
                                FArrayBox&       Utmp,
                                const FArrayBox& a_dV,
                                FArrayBox&  split_tags,
                                BaseFab<unsigned char>& a_Flags,
                                FluxBox&         a_F,
                                Time_Step        *Dts,
                                const Box&       UBox,
                                Grid *grid) 
/*
 *
 *
 *
 *
 ************************************************************* */
{
  CH_assert(isDefined());
  CH_assert(UBox == m_currentBox);

  int nv, ii, jj, kk, nxf, nxb, indf;
  double ***splitcells;
  static Data d;
  static State_1D state;
  static double **u;
  Riemann_Solver *Riemann;
  Index indx;
  Riemann = rsolver;

  nxf = grid[IDIR].np_int + 1;
  nxb = grid[IDIR].lbeg   - 1;

  jj = kk = 0; g_j = jj; g_k = kk;

/* -----------------------------------------------------------------
               Check algorithm compatibilities
   ----------------------------------------------------------------- */

  #if PARABOLIC_FLUX == YES && ARTIFICIAL_VISCOSITY == NO
   print1 ("! Parabolic Terms are not implemented in 1D CTU at the moment\n");
   QUIT_PLUTO(1);
  #endif
 
  #if !(GEOMETRY == CARTESIAN || GEOMETRY == CYLINDRICAL)
   print1 ("! CTU only works in cartesian or cylindrical coordinates\n");
   QUIT_PLUTO(1);
  #endif     

  if (NX1_TOT > NMAX_POINT || NX2_TOT > NMAX_POINT || NX3_TOT > NMAX_POINT){
    print ("! PatchUnsplit: need to re-allocate matrix\n");
    QUIT_PLUTO(1);
  }

/* -----------------------------------------------------------------
                          Allocate memory
   ----------------------------------------------------------------- */

  #if GEOMETRY != CARTESIAN
   for (nv = 0; nv < NVAR; nv++) a_U.divide(a_dV,0,nv);
   #if CHOMBO_CONS_AM == YES
    #if ROTATING_FRAME == YES
     Box curBox = a_U.box();
     for(BoxIterator bit(curBox); bit.ok(); ++bit) {
       const IntVect& iv = bit();
       a_U(iv,iMPHI) /= a_dV(iv,1);
       a_U(iv,iMPHI) -= a_U(iv,RHO)*a_dV(iv,1)*g_OmegaZ;
     }
    #else
     a_U.divide(a_dV,1,iMPHI);
    #endif
   #endif
  #else
   if (g_stretch_fact != 1.) a_U /= g_stretch_fact;
  #endif

  #ifdef SKIP_SPLIT_CELLS
   splitcells = ArrayBoxMap (KBEG, KEND, JBEG, JEND, IBEG, IEND, 
                                 split_tags.dataPtr(0));
  #endif

/* -----------------------------------------------------------
         Allocate static memory areas
   -----------------------------------------------------------  */

  if (state.flux == NULL){
    MakeState (&state);
    d.Vc   = ARRAY_4D(NVAR, 1, 1, NMAX_POINT, double);
    d.flag = ARRAY_3D(1, 1, NMAX_POINT, unsigned char);
    u      = ARRAY_2D(NMAX_POINT, NVAR, double);
  }

  for (ii = 0; ii < NX1_TOT; ii++){
    for (nv = 0; nv < NVAR; nv++){
      u[ii][nv] = a_U.dataPtr(nv)[ii];
    }
    #if COOLING == MINEq || COOLING == H2_COOL
     NormalizeIons(u[ii]);
    #endif
  }

  g_intStage = 1;
  FlagReset (&d);
  getPrimitiveVars (u, &d, grid);
  #if SHOCK_FLATTENING == MULTID
   FindShock (&d, grid);
  #endif
  #ifdef SKIP_SPLIT_CELLS
   for (ii = IBEG; ii <= IEND; ii++){
   if (splitcells[0][0][ii] < 0.5){
     d.flag[0][0][ii] |= FLAG_SPLIT_CELL;
   }}
  #endif

/* ----------------------------------------------------
     Advance solution for a time step g_dt
   ---------------------------------------------------- */

  int numFlux = numFluxes();
  a_F.resize(UBox,numFlux);
  a_F.setVal(0.0);

  g_dir  = IDIR;
  SetIndexes (&indx, grid);
  ResetState (&d, &state, grid);
  for (ii = 0; ii < NX1_TOT; ii++) {
  for (nv = 0; nv < NVAR; nv++) {
    state.v[ii][nv] = d.Vc[nv][0][0][ii];
  }}
  CheckNaN (state.v, IBEG - 1, IEND + 1, 0);
  PrimToCons (state.v, u, 0, NX1_TOT-1);
  States     (&state, IBEG - 1, IEND + 1, grid);
  Riemann   (&state, IBEG - 1, IEND, Dts->cmax, grid);
  RightHandSide (&state, Dts, IBEG, IEND, g_dt, grid);
  saveFluxes (&state, IBEG - 1, IEND, grid);
  for (ii = IBEG; ii <= IEND; ii++) {
  for (nv = 0; nv < NVAR; nv++) {
    u[ii][nv] += state.rhs[ii][nv];
  }}   

// Put fluxes in the FarrayBox a_F to be passed to Chombo

  for (ii = IBEG - 1; ii <= IEND; ii++) {
  for (nv = 0; nv < NVAR; nv++) {
    indf = nv*nxf + ii - nxb;
    a_F[IDIR].dataPtr(0)[indf] = state.flux[ii][nv];
  }}

/* ----------------------------------------------
     Need to include Cooling ? 
     --> convert to primitive space
   ---------------------------------------------- */

  #if COOLING != NO 
   ConsToPrim (u, state.v, IBEG, IEND, state.flag);
   for (ii = IBEG; ii <= IEND; ii++){
   for (nv = NVAR; nv--;  ) {
     d.Vc[nv][0][0][ii] = state.v[ii][nv];
   }} 

/*  ----  Optically thin Cooling sources  ----  */

   #if COOLING == POWER_LAW  /* -- solve exactly -- */
    PowerLawCooling (d.Vc, g_dt, Dts, grid);
   #else
    CoolingSource (&d, g_dt, Dts, grid);
   #endif
 
   for (ii = IBEG; ii <= IEND; ii++){
   for (nv = NVAR; nv--;  ) {
     state.v[ii][nv] = d.Vc[nv][0][0][ii];
   }}

   PrimToCons (state.v, u, IBEG, IEND);
  #endif

  #if CHOMBO_EN_SWITCH == YES
   totEnergySwitch (??, IBEG-IOFFSET, IEND+IOFFSET, 0, 0, 0, 0, -1);
  #endif

/* ------------------------------------
    copy conservative solution vector 
    on to Chombo data structure 
   ------------------------------------ */

/* ---------------------------------------------------------------
    We pass U*dV/m_dx^3 back to Chombo rather than U.
   --------------------------------------------------------------- */

  for (ii = IBEG; ii <= IEND; ii++){
  for (nv = NVAR; nv--;  ) {
    a_U.dataPtr(nv)[ii] = u[ii][nv];
  }} 

  #if GEOMETRY != CARTESIAN
   #if CHOMBO_CONS_AM == YES
    #if ROTATING_FRAME == YES
     for(BoxIterator bit(curBox); bit.ok(); ++bit) {
       const IntVect& iv = bit();
       a_U(iv,iMPHI) += a_U(iv,RHO)*a_dV(iv,1)*g_OmegaZ;
       a_U(iv,iMPHI) *= a_dV(iv,1);
     }
    #else
     a_U.mult(a_dV,1,iMPHI);
    #endif
   #endif
   for (nv = 0; nv < NVAR; nv++) a_U.mult(a_dV,0,nv);
  #else
   if (g_stretch_fact != 1.) a_U *= g_stretch_fact;
  #endif

/* -------------------------------------------------
               Free memory 
   ------------------------------------------------- */

  #ifdef SKIP_SPLIT_CELLS
   FreeArrayBoxMap (splitcells, KBEG, KEND, JBEG, JEND, IBEG, IEND);
  #endif

}

/* ********************************************************* */
void PatchPluto::getPrimitiveVars (double **u, Data *d, Grid *grid)
/*
 *
 *
 * Set physical boundary conditions and convert
 * the conservative matrix U into primitive values
 * d->Vc.
 *
 * 
 *********************************************************** */
{
  int i, j, k, nv, dir;
  int nx;
  int ibeg, iend;
  int lft_side, rgt_side;
  int err;
  static real **v;
  static unsigned char *flag;

  if (v == NULL){
    v    = ARRAY_2D(NMAX_POINT, NVAR, double);
    flag = ARRAY_1D(NMAX_POINT, unsigned char);
  }

  nx = grid[IDIR].np_tot;

/* ------------------------------------------------------- 
     check whether the patch touches a physical boundary
   ------------------------------------------------------- */

  lft_side = (grid[IDIR].lbound != 0);
  rgt_side = (grid[IDIR].rbound != 0);

/* ---------------------------------------------------
     Extract the portion of the domain where U 
     is defined (i.e. NOT in the physical boundary).
   --------------------------------------------------- */

  ibeg = 0; iend = nx - 1;

  #if GEOMETRY == CYLINDRICAL
   for (nv = NVAR; nv--;   ){
   for (i = ibeg; i <= iend; i++){
     u[i][nv] /= grid[IDIR].x[i];
   }}
  #endif

  #if GEOMETRY == SPHERICAL
   for (nv = NVAR; nv--;   ){
   for (i = ibeg; i <= iend; i++){
     u[i][nv] *= m_dx/grid[IDIR].dV[i];
   }}
  #endif

/* -------------------------------------------------
     exclude physical boundaries since the 
     conservative vector U is not yet defined.    
   ------------------------------------------------- */

  if (lft_side) ibeg = IBEG;  
  if (rgt_side) iend = IEND;  

/* ----------------------------------------------
     convert from conservative to primitive 
   ---------------------------------------------- */

  #if CHOMBO_EN_SWITCH == YES
   totEnergySwitch (??, ibeg, iend, 0, 0, 0, 0, +1);
  #endif

  g_dir = IDIR;
  err = ConsToPrim (u, v, ibeg, iend, flag);
  if (err != 0) {
    WARNING(
      print (">> in getPrimitiveVars t = %f\n", g_time);
      showPatch(grid);
      for (i = ibeg; i <= iend; i++){ if (flag[i] != 0) {
        print (">> in zone %d %d\n",i,j);
      }}
    )
  }

  for (i = ibeg; i <= iend; i++){
  for (nv = NVAR; nv--;   ){
    d->Vc[nv][0][0][i] = v[i][nv];
  }}

/* -- now set boundary on primitive variables -- */

  if (lft_side || rgt_side) Boundary (d, ALL_DIR, grid);
}

