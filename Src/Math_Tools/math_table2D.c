/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Miscellaneous functions for handling 2D tables.

  \author A. Mignone (mignone@ph.unito.it)
  \date   June 20, 2014
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

/* ********************************************************************* */
void InitializeTable2D (Table2D *tab, double xmin, double xmax, int nx,
                                      double ymin, double ymax, int ny)
/*!
 * Allocate memory for the arrays contained in the \c *tab structure 
 * and generate a uniformly spaced grid in log10(x), log10(y) within the 
 * range provided by \c xmin, \c xmax, \c ymin, \c ymax.
 * On output, the function initializes the following structure members:
 *
 * - <tt> nx, ny </tt>
 * - <tt> lnxmin, lnxmax, lnymin, lnymax </tt>
 * - <tt> lnx[], lny[] </tt>
 * - <tt> dlnx, dlny  </tt>
 * - <tt> dlnx_1, dlny_1  </tt>
 * - <tt> x[], y[]  </tt>
 * - <tt> dx[], dy[] </tt>
 *
 * \param [in,out] tab   pointer to a Table2D structure
 * \param [in]     xmin  lower column limit. 
 * \param [in]     xmax  upper column limit. 
 * \param [in]     nx    number of equally-spaced column bins (in log space)
 * \param [in]     ymin  lower row limit. 
 * \param [in]     ymax  upper row limit. 
 * \param [in]     ny    number of equally-spaced row bins (in log space)
 * 
 *********************************************************************** */
{
  int i,j;

  tab->nx = nx;
  tab->ny = ny;
    
/* ---------------------------------------------------------------------- 
   If you do not initialize a structure variable, the effect depends on 
   whether it is has static storage (see Storage Class Specifiers) or 
   not. 
   If it is, members with integral types are initialized with 0 and 
   pointer  members are initialized to NULL; otherwise, the value of 
   the structure's members is indeterminate.
  ----------------------------------------------------------------------- */

  tab->f       = ARRAY_2D(tab->ny, tab->nx, double);
  tab->defined = ARRAY_2D(tab->ny, tab->nx, char);
  for (i = 0; i < tab->nx; i++){
    for (j = 0; j < tab->ny; j++) tab->f[j][i] = 0.0;
    for (j = 0; j < tab->ny; j++) tab->defined[j][i] = 1;
  }

  tab->lnx  = ARRAY_1D(tab->nx, double);
  tab->lny  = ARRAY_1D(tab->ny, double);

  tab->x    = ARRAY_1D(tab->nx, double);
  tab->y    = ARRAY_1D(tab->ny, double);

  tab->dx   = ARRAY_1D(tab->nx, double);
  tab->dy   = ARRAY_1D(tab->ny, double);

  tab->dfx  = ARRAY_2D(tab->ny, tab->nx, double);
  tab->dfy  = ARRAY_2D(tab->ny, tab->nx, double);
/*
   -- (beta) function index bookeeping --
    tab->fmin = ARRAY_1D(tab->ny, double);
    tab->fmax = ARRAY_1D(tab->ny, double);
    tab->df   = ARRAY_1D(tab->ny, double);
    tab->nf   = 256;
    tab->i    = ARRAY_2D(tab->ny, tab->nf, int);
*/

/* ---------------------------------------------------------------
    Compute table bounds 
   --------------------------------------------------------------- */

  tab->lnxmin = log10(xmin);
  tab->lnxmax = log10(xmax);
  tab->lnymin = log10(ymin);
  tab->lnymax = log10(ymax);
  
  tab->dlnx = (tab->lnxmax - tab->lnxmin)/(double)(tab->nx - 1.0);
  tab->dlny = (tab->lnymax - tab->lnymin)/(double)(tab->ny - 1.0);

  tab->dlnx_1 = 1.0/tab->dlnx;
  tab->dlny_1 = 1.0/tab->dlny;

  for (i = 0; i < tab->nx; i++) {
    tab->lnx[i] = tab->lnxmin + i*tab->dlnx;
    tab->x[i]   = pow(10.0, tab->lnx[i]);
  }

  for (j = 0; j < tab->ny; j++) {
    tab->lny[j] = tab->lnymin + j*tab->dlny;
    tab->y[j]   = pow(10.0, tab->lny[j]);
  }

/* -- table will not be uniform in x and y: compute spacings -- */

  for (i = 0; i < tab->nx-1; i++) tab->dx[i] = tab->x[i+1] - tab->x[i];
  for (j = 0; j < tab->ny-1; j++) tab->dy[j] = tab->y[j+1] - tab->y[j];

}
/* ********************************************************************* */
void FinalizeTable2D(Table2D *tab)
/*
 *
 *********************************************************************** */
{
  int i,j;
  
/* -------------------------------------------------------
    Compute forward differences 
   ------------------------------------------------------- */
   
  for (j = 0; j < tab->ny; j++){
  for (i = 0; i < tab->nx-1; i++){
    tab->dfx[j][i] = tab->f[j][i+1] - tab->f[j][i];
  }}
  for (j = 0; j < tab->ny-1; j++){
  for (i = 0; i < tab->nx; i++){
    tab->dfy[j][i] = tab->f[j+1][i] - tab->f[j][i];
  }}

}

/* ********************************************************************* */
int LocateIndex(double *yarr, int beg, int end, double y)
/*!
 * Given an array \c yarr[beg..end] and a given value \c y, returns the 
 * array index \c i such that <tt>yarr[i] < y < yarr[i+1]</tt> (for 
 * increasing arrays) or <tt>yarr[i+1] < y < yarr[i]</tt> (for 
 * decreasing arrays).
 * \c yarr must be monotonically increasing or monotonically decreasing.
 *
 * \param [in] yarr   the 1D array 
 * \param [in] beg    initial index of the array
 * \param [in] end    final   index of the array
 * \param [in] y      the value to be searched for
 *
 * \return On success return the index of the array. Otherwise
 *         returns -1.
 *         
 * \b Reference
 *    - "Numerical Recipes in C"
 *
 *********************************************************************** */
{
  int iu,im,il, ascnd;
  
  ascnd = (yarr[end] >= yarr[beg]);  /* ascnd = 1 if table is increasing */

/* --------------------------------------------
    Check if y is bounded between max and min 
   -------------------------------------------- */
   
  if ( (ascnd  && (y < yarr[beg] || y > yarr[end])) ||
       (!ascnd && (y > yarr[beg] || y < yarr[end]))){
/*    print ("! LocatIndex: element outside range ");
    print ("(%12.6e not in [%12.6e, %12.6e]\n",y,yarr[beg],yarr[end]);*/
    return -1;
  }
  
  il = beg;  /*  Initialize lower   */
  iu = end;  /*  and upper limits.  */

  while ( (iu - il) > 1) {     
    im = (iu + il) >> 1;      /* Compute midpoint */
    if (y >= yarr[im] == ascnd) il = im;  /* Replace lower limit   */
    else                        iu = im;  /* or upper limit        */
  }  
  if      (y == yarr[beg]) return beg;
  else if (y == yarr[end]) return end;
  else                     return il;    
}

#define LOG_INTERPOLATION NO
/* ********************************************************************* */
int InverseLookupTable2D (Table2D *tab, double y, double f, double *x)
/*!
 * Perform inverse lookup table interpolation: 
 * given a 2D table <tt>f(i,j)=f(x(i), y(j))</tt>
 * and the value of \c y, find the value of \c x.
 * The algorithm proceeds by a combination of 1D lookup table algorithms:
 *
 * - Find the value of \c j such that <tt>y[j]<y<y[j+1]</tt>;
 * - The solution must then lie on the intermediate line between
 *   \c j and \c j+1 where the value of \c y is given.
 * - To find the horizontal coordinate, we first locate initial and final 
 *   indices \c ib and \c ie at \c j and \c j+1 using 1D lookup table.
 * - Then construct a 1D array between these two indices values, perform
 *   again 1D lookup table to find the actual index \c i and invert 
 *   the bilinear interpolant by solving
 *   \f[
 *      f =   f_{i,j}(1 - x_n)(1 - y_n) + f_{i+1,j}x_n(1 - y_n)
 *          + f_{i,j+1}(1 - x_n)y_n     + f_{i+1,j+1}x_ny_n;
 *   \f]
 *   for \f$ x_n\f$ (normlized coordinate between 0 and 1).
 *
 * \note The table must be monotonic in both \c and \c y.
 *
 * \param [in]   tab   a pointer to a Table2D structure
 * \param [in]   y     the ordinata 
 * \param [in]   f     the value of the function
 * \param [out]  *x    the value of the abscissa.
 *
 * \return  Returns 0 on success, otherwise:
 *          - -2 if \c y is below range; 
 *          -  2 if \c y is above range.
 *          
 *********************************************************************** */
{
  int    i, j, k, i0, i1, ib, ie;
  double lny, xn, yn;
  double **ftab, **dfx, **dfy;
  static double *f1;
 
  if (f1 == NULL) f1 = ARRAY_1D(8192, double);

  ftab = tab->f;
  dfx  = tab->dfx;
  dfy  = tab->dfy;

  lny = log10(y);
  if (lny < tab->lnymin){
    print ("! InverseLookupTable2D: lny outside range: %12.6e < %12.6e\n",
            lny, tab->lnymin);
    return -2;
  }
  if (lny > tab->lnymax){
    print ("! InverseLookupTable2D: lny outside range: %12.6e > %12.6e\n",
            lny, tab->lnymax);
    return 2;
  }

  j   = INT_FLOOR((lny - tab->lnymin)*tab->dlny_1); /* Find row */
  #if LOG_INTERPOLATION == YES
   yn  = (lny - tab->lny[j])/tab->dlny;
  #else
   yn  = (y - tab->y[j])/tab->dy[j];
  #endif

  i0 = LocateIndex(tab->f[j],   0, tab->nx-1, f);
  i1 = LocateIndex(tab->f[j+1], 0, tab->nx-1, f);

  if (i0 < 0 || i1 < 0) return 2;

  ib = MIN(i0,i1);
  ie = MAX(i0,i1) + 1;

  if (ie > (ib+1)){  
    for (i = ib; i <= ie; i++) {
      f1[i] = tab->f[j][i]*(1.0 - yn) + tab->f[j+1][i]*yn;
    }
    i = LocateIndex(f1, ib, ie, f);
  }else{
    i = ib; 
  }

  xn  = f - ftab[j][i] - yn*dfy[j][i];
  xn /= yn*(dfx[j+1][i] - dfx[j][i]) + dfx[j][i];

  #if LOG_INTERPOLATION == YES
   (*x) = xn*tab->dlnx + tab->lnx[i];
   (*x) = pow(10.0, *x);
  #else
   (*x) = xn*tab->dx[i] + tab->x[i];
  #endif

/* -- use for verification 
double fint;
Table2DInterpolate(tab, *x, y, &fint);
if (fabs(f/fint-1.0) > 1.e-6){
  printf ("! InverseLookupTable2D: direct table does not match\n");
  QUIT_PLUTO(1);
}

if ( y < tab->y[j] || y > tab->y[j+1]){
  printf ("!err in InverseLookupTable2D\n");
  QUIT_PLUTO(1);
}
if ( f > ftab[j][i] && f < ftab[j+1][i+1]) return 0;
else {
  printf ("!err in InverseLookupTable2D\n");
  QUIT_PLUTO(1);
}
*/
  return 0;

/*
if (   tab->defined[j][i] && tab->defined[j+1][i]
    && tab->defined[j][i+1] && tab->defined[j+1][i+1]){
  (*x) = ...
  return 0;   
}else return 1;
*/
}  

/* ********************************************************************* */
int Table2DInterpolate(Table2D *tab, double x, double y, double *f)
/*!
 * Use bilinear interpolation to find the function f(x,y) from the 
 * 2D table \c *tab.
 * Since the grid is equally spaced in \c log(x) and \c log(y), 
 * the (i,j) indices such that 
 * <tt> x[i] <= x < x[i+1] </tt> and <tt> y[j] <= y < y[j+1] </tt>
 * are found by a simple division.
 * Then bilinear interpolation can be done in either log or linear
 * coordinates. The latter is expected to be slightly faster since it
 * avoids one \c pow() operation.
 * 
 * \param [in]   *tab   a pointer to a Table2D structure
 * \param [in]   x      the abscissa where interpolation is needed.
 * \param [in]   y      the ordinata where interpolation is needed.
 * \param [out] *f      the interpolated value
 *
 * \return  - Return 0 on success, otherwise: 
 *          - -1 if \c x is below column range; 
 *          -  1 if \c x is above column range;
 *          - -2 if \c y is below row range; 
 *          -  2 if \c y is above row range.
 *
 *********************************************************************** */
{
  int i,j;
  double lnx, lny, xn, yn;
  
  lnx = log10(x);  /* Take the log to find the indices */
  lny = log10(y);

/* -------------------------------------------------------
    Check bounds 
   ------------------------------------------------------- */

  if (lnx < tab->lnxmin){
    WARNING(
      print ("! Table2DInterpolate: lnx outside range %12.6e < %12.6e\n",
              lnx, tab->lnxmin);
    )
    return -1;
  }
   
  if (lnx > tab->lnxmax){
    WARNING(
      print ("! Table2DInterpolate: lnx outside range: %12.6e > %12.6e\n",
              lnx, tab->lnxmax);
    )
    return 1;
  }

  if (lny < tab->lnymin){
    WARNING(
      print ("! Table2DInterpolate: lny outside range: %12.6e < %12.6e\n",
              lny, tab->lnymin);
    )
    return -2;
  }

  if (lny > tab->lnymax){
    WARNING(
      print ("! Table2DInterpolate: lny outside range: %12.6e > %12.6e\n",
              lny, tab->lnymax);
    )
    return 2;
  }

/* ------------------------------------------------
    Find column and row indices i and j 
   ------------------------------------------------ */
   
  i  = INT_FLOOR((lnx - tab->lnxmin)*tab->dlnx_1); 
  j  = INT_FLOOR((lny - tab->lnymin)*tab->dlny_1); 

  #if LOG_INTERPOLATION == YES
   xn = (lnx - tab->lnx[i])/tab->dlnx;  /* Compute normalized log coordinates */ 
   yn = (lny - tab->lny[j])/tab->dlny;  /* on the unit square [0,1]x[0,1]    */
  #else
   xn = (x - tab->x[i])/tab->dx[i];  /* Compute normalized linear coordinates */ 
   yn = (y - tab->y[j])/tab->dy[j];  /* on the unit square [0,1]x[0,1]        */
  #endif
  
  (*f) =    tab->f[j][i]*(1.0 - xn)*(1.0 - yn)
          + tab->f[j][i+1]*xn*(1.0 - yn)
          + tab->f[j+1][i]*(1.0 - xn)*yn
          + tab->f[j+1][i+1]*xn*yn;
  return 0;  /* success */
}

/* ********************************************************************* */
void WriteBinaryTable2D (char *fname, Table2D *tab)
/*! 
 *  The binary table is a compact format used to write a 2D array 
 *  together with simple structured coordinates.
 *  The file consists of the following information:
 *  \verbatim
     nx
     ny
     <x[0]..x[nx-1]>
     <y[0]..y[ny-1]>
     <q[0][0]..q[0][nx-1]
      q[1][0]..q[1][nx-1]
        .......
      q[ny-1][0]..q[ny-1][nx-1]>
    \endverbatim
 *  All fields are written in binary format using double precision
 *  arithmetic.
 *
 *********************************************************************** */
{
  int    i,j;
  double scrh;
  FILE *fp;
  
  fp = fopen(fname,"wb");
 
  scrh = (double)tab->nx;
  fwrite(&scrh, sizeof(double), 1, fp);
  
  scrh = (double)tab->ny;
  fwrite(&scrh, sizeof(double), 1, fp);
  
  fwrite(tab->lnx, sizeof(double), tab->nx, fp);
  fwrite(tab->lny, sizeof(double), tab->ny, fp);
  
  for (j = 0; j < tab->ny; j++){
    fwrite (tab->f[j], sizeof(double), tab->nx, fp);
  }
  
  fprintf (fp,"\n");
  fclose(fp);
}
