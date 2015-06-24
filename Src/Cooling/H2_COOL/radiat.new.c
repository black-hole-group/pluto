#include "pluto.h"

#define POW10(x)  (pow(10.0, x))
#define TABLE_NPT   2048
static void Make_RateTables(double, double *);

/* ********************************************************************* */
void Make_RateTables(double T, double *krvals)
/*
 *********************************************************************** */
{
  int i, indx_lo, indx_hi;
  double st, t3, tev, lnT, lnTmin, lnTmax, dlnT;
  double scrh, dum_lo, dum_hi;
  double Tmin, Tmax;
  static double *lnTarr, *Rate_arr1, *Rate_arr2, *Rate_arr3, *Rate_arr4;
  static double *Rate_arr5, *Rate_arr6;
  static double *Emiss_arr1, *Emiss_arr2, *Emiss_arr3;
  double Tdust = 15.0;

  Tmin = 1.0e-2;
  Tmax = 1.0e8;
  lnTmin = log10(Tmin);
  lnTmax = log10(Tmax);
  dlnT   = (lnTmax - lnTmin)/((double)TABLE_NPT - 1.0);

  if (T < Tmin || T > Tmax){
    print1 ("!Radiat: Temperature %12.6e for Rate Tables out of range.\n", T);
    QUIT_PLUTO(1);
  }

  if (lnTarr == NULL){
    lnTarr     = ARRAY_1D(TABLE_NPT, double);
    Rate_arr1  = ARRAY_1D(TABLE_NPT, double);
    Rate_arr2  = ARRAY_1D(TABLE_NPT, double);
    Rate_arr3  = ARRAY_1D(TABLE_NPT, double);
    Rate_arr4  = ARRAY_1D(TABLE_NPT, double);
    Rate_arr5  = ARRAY_1D(TABLE_NPT, double);
    Rate_arr6  = ARRAY_1D(TABLE_NPT, double);
    Emiss_arr1 = ARRAY_1D(TABLE_NPT, double);
    Emiss_arr2 = ARRAY_1D(TABLE_NPT, double);
    Emiss_arr3 = ARRAY_1D(TABLE_NPT, double);
    
    for(i = 0 ; i < TABLE_NPT; i++){
      lnTarr[i] = lnTmin + i*dlnT;
      T     = POW10(lnTarr[i]);
      st    = sqrt(T);
      t3    = T/1.e3;  
      tev   = T*CONST_kB/CONST_eV;
      
      Rate_arr1[i] = 3.e-18*st/(1.0 + 0.04*st + 2.e-3*T + 8.e-6*T*T); /* fa ~ 1 and Tg << T  */
      Rate_arr2[i] = 1.067e-10*pow(tev,2.012)*exp(-(4.463/tev)*pow((1.0 + 0.2472),3.512));
      Rate_arr3[i] = 1.e-8*exp(-84100.0/T); 
      Rate_arr4[i] = 4.4e-10*pow(T,0.35)*exp(-102000.0/T);
      Rate_arr5[i] = 2.6e-11/st;
      Rate_arr6[i] = 1.08e-8*st*exp(-157890.0/T)/(13.6*13.6);
      Emiss_arr1[i] = 6.7e-19*exp(-5.86/t3) + 1.6e-18*exp(-11.7/t3);
      Emiss_arr2[i] = (9.5e-22*pow(t3,3.76)/(1. + 0.12*pow(t3,2.1)))*exp(pow((-0.13/t3),3.0))
	+ 3.0e-24*exp(-0.51/t3); 
      Emiss_arr3[i] = 3.8e-33*st*(T-Tdust)*(1.0 - 0.8*exp(-75.0/T));
    }
  }else{   /* Perform linear interpolation */
    lnT = log10(T);

    scrh    = (lnT - lnTmin)/dlnT;
    indx_lo = INT_FLOOR(scrh);
    indx_hi = indx_lo + 1;
if (lnT > lnTarr[indx_hi] || lnT < lnTarr[indx_lo]){
  print ("!!!!!!! this should never happen  %12.6e\n", scrh);
  QUIT_PLUTO(1);
}
    dum_lo  = (lnTarr[indx_hi] - lnT)/dlnT;
    dum_hi  = (lnT - lnTarr[indx_lo])/dlnT;

    if (indx_lo < 0){
      print ("! Make_RateTable: Cannot Locate Index for T = %12.6e\n",T);
      QUIT_PLUTO(1);
    }else{
      krvals[0] = Rate_arr1[indx_lo]*dum_lo + Rate_arr1[indx_hi]*dum_hi; /* kr1 */
      krvals[1] = Rate_arr2[indx_lo]*dum_lo + Rate_arr2[indx_hi]*dum_hi; /* kr2 */
      krvals[2] = Rate_arr3[indx_lo]*dum_lo + Rate_arr3[indx_hi]*dum_hi; /* kr3 */
      krvals[3] = Rate_arr4[indx_lo]*dum_lo + Rate_arr4[indx_hi]*dum_hi; /* kr4 */
      krvals[4] = Rate_arr5[indx_lo]*dum_lo + Rate_arr5[indx_hi]*dum_hi; /* cr */
      krvals[5] = Rate_arr6[indx_lo]*dum_lo + Rate_arr6[indx_hi]*dum_hi; /* ci */
      krvals[6] = Emiss_arr1[indx_lo]*dum_lo + Emiss_arr1[indx_hi]*dum_hi; /* em[1] */
      krvals[7] = Emiss_arr2[indx_lo]*dum_lo + Emiss_arr2[indx_hi]*dum_hi; /* em[2] */
      krvals[8] = Emiss_arr3[indx_lo]*dum_lo + Emiss_arr3[indx_hi]*dum_hi; /* em[11] */
    }
  }
  return;
}

/* ***************************************************************** */
void Radiat (double *v, double *rhs)
/*
 *  Cooling for neutral or singly ionized gas: good up to about 35,000 K
 *  in equilibrium or shocks in neutral gas up to about 80 km/s.
 *  Assumed abundances in ab
 *  Uses t : Kelvin
 *       dene :  electron density cm*-3
 *       fneut : hydrogen neutral fraction    (adimensionale)
 *       ci,cr : H ionization and recombination rate coefficients 
 *
 *
 * em(1) = TOTAL EMISSIVITY : (ergs cm**3 s**-1)
 * em(2) = Ly alpha  + two photon continuum: Aggarwal MNRAS 202,
 *         10**4.3 K
 * em(3) = H alpha: Aggarwal, Case B
 * em(4) = He I 584 + two photon + 623 (all n=2 excitations): Berrington
 *         &Kingston,JPB 20
 * em(5) = C I 9850 + 9823: Mendoza, IAU 103, 5000 K
 * em(6) = C II, 156 micron: Mendoza, 10,000 K
 * em(7) = C II] 2325 A: Mendoza, 15,000 K
 * em(8) = N I 5200 A: Mendoza, 7500 K
 * em(9) = N II 6584 + 6548 A: Mendoza
 * em(10) = O I 63 micron: Mendoza,2500 K
 * em(11) = O I 6300 A + 6363 A: Mendoza, 7500 K
 * em(12) = O II 3727: Mendoza
 * em(13) = Mg II 2800: Mendoza
 * em(14) = Si II 35 micron: Dufton&Kingston, MNRAS 248
 * em(15) = S II 6717+6727: Mendoza
 * em(16) = Fe II 25 micron: Nussbaumer&Storey
 * em(17) = Fe II 1.6 micron
 * em(18) = thermal energy lost by ionization
 * em(19) = thermal energy lost by recombination (2/3 kT per
 *          recombination.
 *          The ionization energy lost is not included here.
 *
 ******************************************************************* */
{
  int   ii, k, nv, status;
  double  T, mu, t3, rho, prs, fn, gn, hn, x;
  double  N_H, n_el, cr, ci, rlosst, src_pr, crold,crnew;
  double kr1, kr2, kr3, kr4, em[20], krvalues[9];
  
  static int first_call = 1;
  static real E_cost, Unit_Time, N_H_rho, frac_He, frac_Z, N_tot;

  if (first_call) {

    E_cost    = UNIT_LENGTH/UNIT_DENSITY/pow(UNIT_VELOCITY, 3.0);
    Unit_Time = UNIT_LENGTH/UNIT_VELOCITY;
        
    N_H_rho  = (UNIT_DENSITY/CONST_amu)*(H_MASS_FRAC/CONST_AH);
    frac_He  = (He_MASS_FRAC/CONST_AHe)*(CONST_AH/H_MASS_FRAC);
    frac_Z   = ((1 - H_MASS_FRAC - He_MASS_FRAC)/CONST_AZ)*(CONST_AH/H_MASS_FRAC);
    N_tot  =  N_H_rho*(1.0 + frac_He + frac_Z);
    Make_RateTables(100.0, krvalues);
    first_call = 0;
  }

/* ---------------------------------------------
    Force fneut and fmol to stay between [0,1]
   --------------------------------------------- */

  for (nv = NFLX; nv < (NFLX + NIONS); nv++){
    v[nv] = MAX(v[nv], 0.0);
    v[nv] = MIN(v[nv], 1.0);
  }

  v[X_H2] = MIN(v[X_H2], 0.5);

/* ---------------------------------------------------
    Compute temperature from rho, rhoe and fractions
   --------------------------------------------------- */
   
  rho = v[RHO];
  mu  = MeanMolecularWeight(v); 
  if (mu < 0.0){
    print1 ("! Radiat: mu = %f < 0 \n",mu);
    QUIT_PLUTO(1);
  }
  #if EOS == IDEAL
   if (v[RHOE] < 0.0) v[RHOE] = g_smallPressure/(g_gamma-1.0);
   prs = v[RHOE]*(g_gamma-1.0);
   T   = prs/rho*KELVIN*mu;
   T   = MAX(T, g_minCoolingTemp);
  #else
   status = GetEV_Temperature(v[RHOE], v, &T);
   if (status != 0) {
     T = T_CUT_RHOE;
     v[RHOE] = InternalEnergy(v, T);
   }
  #endif

  N_H = N_H_rho*rho;  /* Hydrogen number density
                         N_H = N(X_HI) + 2N(X_H2)+  N(X_HII)  */
  fn  = v[X_HI];
  gn  = v[X_H2];
  hn  = v[X_HII];

/* Recombination and ionization coefficients */ 
  

  n_el = N_H*(hn + 0.5*CONST_AZ*frac_Z);  /* -- electron number density, in cm^{-3} -- */
  
  if (n_el < 0){
    print1 ("!NEG ELECTRONS --No of Electrons is negative");
    QUIT_PLUTO(1);
  }

  t3    = T/1.e3;  
 
/* The Units of kr1, kr2, kr3, kr4, cr, ci = cm^{3}s^{-1}*/
  
  Make_RateTables(T, krvalues);
  kr1 = krvalues[0];
  kr2 = krvalues[1];
  kr3 = krvalues[2];
  kr4 = krvalues[3];
  cr  = krvalues[4];
  ci  = krvalues[5];

  /*
  kr1 = 3.e-18*st/(1.0 + 0.04*st + 2.e-3*T + 8.e-6*T*T);
  kr2 = 1.067e-10*pow(tev,2.012)*exp(-(4.463/tev)*pow((1.0 + 0.2472),3.512));
  kr3 = 1.e-8*exp(-84100.0/T); 
  kr4 = 4.4e-10*pow(T,0.35)*exp(-102000.0/T);
 */
/* ************ ENSURING SUM IS UNITY ********************** 
  * This is done to ensure that the sum of speicies 
  * equal to 1. Total Hydrogen number N_H composed of
  * N_H = nHI + 2.*nH2 + nHII. 
  * while fraction gn = nH2/N_H, thus the rhs[X_H2] has 
  * a pre-factor of 0.5. 
  * And the fact the sum = 1 is ensured if the
  * corresponding sum of RHS = 0. 
  * i.e., rhs[X_H1] + 2.0*rhs[X_H2] + rhs[X_HII] = 0.0.
  ********************************************************** */

  x = n_el/N_H;
  rhs[X_HI] = Unit_Time*N_H*(  gn*(kr2*fn + kr3*gn + kr4*x) + cr*hn*x 
                              -fn*(ci*x + kr1*fn));
  rhs[X_H2]  = 0.5*Unit_Time*N_H*(kr1*fn*fn - gn*(kr2*fn + kr3*gn + kr4*x));
  rhs[X_HII] = Unit_Time*N_H*(ci*fn - cr*hn)*x;
/*
  double sum;
  sum = 0.0;
  sum=fabs(rhs[X_HI] + 2.0*rhs[X_H2] + rhs[X_HII]);
  if (sum > 1.e-9){
    print("%le > Sum(RHS) != 0 for H!!\n",sum);
    QUIT_PLUTO(1);
  }
*/

  double logt3 = log(t3);

/* The unit of cooling function Lambda in erg cm^{-3} s^{-1}*/
  
  /* em[1] = 6.7e-19*exp(-5.86/t3) + 1.6e-18*exp(-11.7/t3);
  em[2] = (9.5e-22*pow(t3,3.76)/(1. + 0.12*pow(t3,2.1)))*exp(pow((-0.13/t3),3.0))
           + 3.0e-24*exp(-0.51/t3);*/
  

  em[1] = krvalues[6];
  em[2] = krvalues[7];
  em[3] = em[1] + em[2];

  if (t3 < 0.2){ /* Following the Table 8 from Glover & Abel, MNRAS (2008) */
    em[4] = - 16.818342 + 399.521483*logt3; 
//+ 37.383713*logt3 + 58.145166*2.0*logt3 + 48.656103*3.0*logt3 + 20.159831*4.0*logt3 + 3.8479610*5.0*logt3;
  }else if(t3 < 1.0){
    em[4] = - 24.311209 - 209.2192867*logt3;
//+ 3.5692468*logt3 - 11.332860*2.0*logt3 - 27.850082*3.0*logt3 - 21.328264*4.0*logt3 - 4.2519023*5.0*logt3;
  }else{
    em[4] = - 24.311209 + 0.7397324*logt3;
//+ 4.6450521*logt3 - 3.7209846*2.0*logt3 + 5.9369081*3.0*logt3 - 5.5108047*4.0*logt3 + 1.5538288*5.0*logt3;
  }
  em[4] = MIN(em[4],50.0);  /* Avoid overflow problems when taking pow() */
  em[4] = POW10(em[4]);

  em[5] = - 23.962112 + 1.097389*logt3;
//+ 2.09433740*logt3 - 0.77151436*2.0*logt3 + 0.43693353*3.0*logt3 - 0.14913216*4.0*logt3 - 0.033638326*5.0*logt3; 

  em[5] = POW10(em[5]);
  em[6] = (fn*em[4] + gn*em[5])*N_H; 
  
/* Cooling due to X_H2 rotational vibration.  */

  em[7] = gn*N_H*em[3]/(1.0 + (em[3]/em[6])); 

/* Cooling due to ionization */
  
  em[8] = ci*13.6*1.6e-12*fn;           
  
/* Cooling due to radiative recombination */

  em[9] = cr*0.67*1.6e-12*hn*T/11590.0; 
  
/* Cooling due to H2 disassociation */

  em[10] = 7.18e-12*(kr2*fn*gn + kr3*gn*gn + kr4*gn*(n_el/N_H))*N_H*N_H; 

/* Cooling due to gas grain cooling */

  em[11] = krvalues[8]*N_H*N_H; //3.8e-33*st*(T-Tdust)*(1.0 - 0.8*exp(-75.0/T))*N_H*N_H;

  em[13]  = em[11] + em[10] + em[7] + (em[8] + em[9])*n_el*N_H;    
     
/* ---------------------------------------------------
    rlosst is the energy loss in units of erg/cm^3/s;
    it must be multiplied by cost_E in order to match 
    non-dimensional units.
    Source term for the neutral fraction scales with 
    UNIT_TIME
   --------------------------------------------------- */

  rlosst  =  em[13];
  rhs[RHOE] = -E_cost*rlosst;
  rhs[RHOE] *= 1.0/(1.0 + exp(-(T - g_minCoolingTemp)/100.0)); /* -- lower cutoff -- */
}

/* ********************************************************************* */
double MeanMolecularWeight (double *V)
/*
 *   Compute the mean molecular weight as function of the 
 *   composition of the gas.
 *   The definitiion of the mean molecular weight \mu is 
 *   the standard one:
 *
 *     1     \sum_k f_k n_k
 *    --- = ----------------     (Clayton, pag 82-83)
 *    \mu    \sum_k f_k A_k
 * 
 *   where 
 *
 *    f_k   : is the fractional abundance (by number) with
 *            respect to hydrogen, f_k = N_k/N_H
 *
 *    A_K   : is the atomic weight
 *
 *    n_k   : is the number of free particles 
 *            contributed to the gas by element k
 *
 *   The mean molecular weight satifies 
 *
 *               \rho = \mu m_{amu} N_{tot}
 *   
 *   where N_{tot} is the total number of particles
 *
 *   For the ``Raymond'' cooling module \mu is calculated
 *   as follows:
 *
 *            A_H + f_He*A_He + f_Z*A_z
 *    \mu =  ---------------------------
 *             2 - fn + f_He + 2*f_Z
 *
 * 
 * ARGUMENTS
 *
 *   V:   a set of primitive variables
 *
 *********************************************************************** */
{
  int nv;
  double rho_mH, munum, muden;
 
  for (nv = NFLX; nv < (NFLX + NIONS); nv++){
    V[nv] = MAX(V[nv], 0.0);
    V[nv] = MIN(V[nv], 1.0);
  }

  V[X_H2] = MIN(V[X_H2], 0.5);
  
  rho_mH = V[RHO]*(UNIT_DENSITY/CONST_amu);
  double N_H  = rho_mH*(H_MASS_FRAC/CONST_AH);  
  double N_He = rho_mH*(He_MASS_FRAC/CONST_AHe); 
  double N_Z  = rho_mH*((1.0-H_MASS_FRAC-He_MASS_FRAC)/CONST_AZ);  

  double fracHe = N_He/N_H;
  double fracZ = N_Z/N_H;
  
  double fn = V[X_HI];
  double gn = V[X_H2];
  double hn = V[X_HII];
 
  munum = 1.0 + CONST_AHe*(fracHe) + CONST_AZ*(fracZ);
  muden = fn + gn + 2*hn + fracHe + fracZ + 0.5*CONST_AZ*(fracZ);

  return munum/muden;
}
