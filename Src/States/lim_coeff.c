#include"pluto.h"

static double **ap, **am;

void MAKE_LIM_COEFF (Grid *grid);
{

  if (ap == NULL){ 

    ap = ARRAY_2D(3, NMAX_POINT, NVAR);
    am = ARRAY_2D(3, NMAX_POINT, NVAR);

    dx = grid[IDIR].dx;
    dy = grid[JDIR].dx;
    dz = grid[KDIR].dx;

    for (i = 1; i < NX1_TOT-1; i++){
     ap[IDIR][i] = dx[i]/(dx[i] + dx[i + 1]);
     am[IDIR][i] = dx[i]/(dx[i] + dx[i - 1]);
    }
    for (i = 1; i < NX2_TOT-1; i++){
     ap[JDIR][i] = dy[i]/(dy[i] + dy[i + 1]);
     am[JDIR][i] = dy[i]/(dy[i] + dy[i - 1]);
    }
    for (i = 1; i < NX3_TOT-1; i++){
     ap[KDIR][i] = dz[i]/(dz[i] + dz[i + 1]);
     am[KDIR][i] = dz[i]/(dz[i] + dz[i - 1]);
    }
  }
   

}


double *GET_AP_COEFF(int dir)
{
  return (ap[dir]);
}

