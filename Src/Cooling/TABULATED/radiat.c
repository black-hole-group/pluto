#include "pluto.h"
/* ***************************************************************** */
void Radiat (double *v, double *rhs)
/*!
 *   Provide r.h.s. for tabulated cooling.
 * 
 ******************************************************************* */

{
	
	int klo, khi, kmid, repvar, index;
	static int length_cool, length_n, length_radius, length_T;
	double dn, nlo, nhi, nmid, n;
	double dradius, radiuslo, radiushi, radiusmid, radius;
	double dT, Tlo, Thi, Tmid, scrh, T;
	double c1, c2, c3, c4, c5, c6, c7, prs, mu_e, k1, mu;
	
	static double *T_tab, *n_tab, *radius_tab, *cool_tab, E_cost;

	FILE *fcool, *ftemp, *fndens, *fradius;

	if (T_tab == NULL || n_tab == NULL || radius_tab == NULL){
    		print1 (" > Reading table from disk...\n");
    		fcool = fopen("coolingtable.dat","r");
		ftemp = fopen("temperature.dat","r");
		fndens = fopen("eletronicdensity.dat","r");
		fradius = fopen("radius.dat","r");

    		if (fcool == NULL || ftemp == NULL || fndens == NULL || fradius == NULL){
      			print1 ("! Radiat: one or more of the tables could not be found.\n");
     			QUIT_PLUTO(1);
    		}
    		cool_tab = ARRAY_1D(1000000, double);
    		T_tab = ARRAY_1D(100, double);
		n_tab = ARRAY_1D(100, double);
		radius_tab = ARRAY_1D(100, double);

    		length_T = 0;
    		while (fscanf(ftemp, "%lf\n", T_tab + length_T)!=EOF) {
     	 		length_T++;
   		}

		length_n = 0;
    		while (fscanf(fndens, "%lf\n", n_tab + length_n)!=EOF) {
     	 		length_n++;
   		}

		length_radius = 0;
    		while (fscanf(fradius, "%lf\n", radius_tab + length_radius)!=EOF) {
     	 		length_radius++;
   		}

		length_cool = 0;
    		while (fscanf(fcool, "%lf\n", cool_tab + length_cool)!=EOF) {
     	 		length_cool++;
   		}
   		E_cost = UNIT_LENGTH/UNIT_DENSITY/pow(UNIT_VELOCITY, 3.0); 	/*Conversion constant of cooling rate (cm3 s / erg)*/
  	}


	/* ---------------------------------------------
       	    Get pressure, temperature, number density
	    and distance	 
   	--------------------------------------------- */


	prs = v[RHOE]*(g_gamma-1.0);						/*Eq 7.6 PLUTO*/
  	if (prs < 0.0) {							/*If pressure is < 0 we set a minimum pressure to calculate v[RHOE]*/
    		prs     = g_smallPressure;
    		v[RHOE] = prs/(g_gamma - 1.0);
  	}

  	mu  = MeanMolecularWeight(v);
  	T   = prs/v[RHO]*KELVIN*mu;						/*Eq 7.2 PLUTO*/

	if (T != T){								/*Verify if temperature is NaN*/
    		printf (" ! Nan found in radiat \n");
    		printf (" ! rho = %12.6e, prs = %12.6e\n",v[RHO], prs);
    		QUIT_PLUTO(1);
  	}

	/* ---------------------------------------------
       	    	    Get electron temperature	 
   	--------------------------------------------- */
	
	k1 = 1.0;
	T = T / (pow(v[VX1] / 100.0, -k1) + 2.0);				/*Electron Temperature according to Wu et al. 2016*/

  	if (T < g_minCoolingTemp) { 						/*Verify if it has reached the minimum temperature*/
    		rhs[RHOE] = 0.0;
    		return;
  	}
	
	radius = v[VX1] * UNIT_LENGTH;
	mu_e = 2/(1 + H_MASS_FRAC);						/*Effective mean molecular weight of electrons for fully ionized gas*/
	n = v[RHO] * UNIT_DENSITY / (mu_e * CONST_amu);


	/* ----------------------------------------------
        	Table lookup by binary search  
   	---------------------------------------------- */


	klo = 0;
	khi = length_n - 1;	

	while (klo != (khi - 1)) {						
   		kmid = (klo + khi) / 2;						
   		nmid = n_tab[kmid];						
    		if (n <= nmid) {						
     		khi = kmid;							
    		}else if (n > nmid) {						
 			klo = kmid;
    		 }
  	}

	dn = n_tab[khi] - n_tab[klo];
	nlo = n_tab[klo];
	nhi = n_tab[khi];

	repvar = length_cool / length_n;
	index = repvar * klo;

	klo = 0;
	khi = length_radius - 1;
	
	while (klo != (khi - 1)) {						
   		 kmid = (klo + khi) / 2;						
    		 radiusmid = radius_tab[kmid];						
    		 if (radius <= radiusmid) {						
     		 khi = kmid;							
    		 }else if (radius > radiusmid) {						
 		  	klo = kmid;
    		  }
  	}
	
	dradius = radius_tab[khi] - radius_tab[klo];
	radiuslo = radius_tab[klo];
	radiushi = radius_tab[khi];

	repvar = repvar / length_radius;
	index += repvar * klo;

	klo = 0;
	khi = length_T - 1;

	while (klo != (khi - 1)) {						
   		 kmid = (klo + khi) / 2;						
    		 Tmid = T_tab[kmid];						
    		 if (T <= Tmid) {						
     		 khi = kmid;							
    		 }else if (T > Tmid) {						
 		  	klo = kmid;
    		  }
  	}

	dT = T_tab[khi] - T_tab[klo];
	Tlo = T_tab[klo];
	Thi = T_tab[khi];
 	
	index += klo;


	/*---------------------------------
	   Weighing the cooling rate value
	  ---------------------------------*/
	

	c1 = dT * dn * dradius;
	c2 = Tlo - T;
	c3 = radiuslo - radius;
	c4 = nlo - n;
	c5 = T - Thi;
	c6 = radius - radiushi;
	c7 = n - nhi;

	scrh = cool_tab[index] * c2 * c3 * c4 / c1;
	index++;
	scrh += cool_tab[index] * c5 * c3 * c4 / c1;
	index += - 1 + length_T;
	scrh += cool_tab[index] * c2 * c6 * c4 / c1;
	index++;
	scrh += cool_tab[index] * c5 * c6 * c4 / c1;
	index += length_T * (length_radius - 1) - 1;
	scrh += cool_tab[index] * c2 * c3 * c7 / c1;
	index++;
	scrh += cool_tab[index] * c5 * c3 * c7 / c1;
	index += length_T - 1;
	scrh += cool_tab[index] * c2 * c6 * c7 / c1;
	index++;
	scrh += cool_tab[index] * c5 * c6 * c7 / c1;


	/* -----------------------------------------------
    		Compute r.h.s
   	----------------------------------------------- */

											
  	rhs[RHOE] = -scrh*E_cost;							/*Eq 9.1 and converting cooling rate to code units*/

}
