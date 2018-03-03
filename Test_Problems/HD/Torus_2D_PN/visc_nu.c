/* /////////////////////////////////////////////////////////////////// */
/*! \file
 *  \brief Specification of explicit first and second viscosity coefficients*/
/* /////////////////////////////////////////////////////////////////// */
#include "pluto.h"
/* ************************************************************************** */
void Visc_nu(double *v, double x1, double x2, double x3,
                        double *nu1, double *nu2)
/*!
 *
 *  \param [in]      v  pointer to data array containing cell-centered quantities
 *  \param [in]      x1 real, coordinate value
 *  \param [in]      x2 real, coordinate value
 *  \param [in]      x3 real, coordinate value
 *  \param [in, out] nu1  pointer to first viscous coefficient
 *  \param [in, out] nu2  pointer to second viscous coefficient
 *
 *  \return This function has no return value.
 * ************************************************************************** */
{
/*
Here we adopt the alpha prescription to mimic the MRI-driven viscosity.
Our implementation follows Stone et al. (1999); Yuan et al. (2012). 

Stone+1999 proposed two different prescriptions for the kinematic 
viscosity (second paragraph of page 1004):

1) nu = alpha * rho => C: *nu1 = alpha * v[RHO] * v[RHO]
2) nu = alpha * sqrt(r) => C: *nu1 = alpha * v[RHO] * sqrt(x1)

To allow a more direct comparison with Yuan+2012, we will adopt the
second viscosity prescription.
*/
  double alpha = 0.3; // Shakura-Sunyaev alpha parameter
  *nu1 = alpha * v[RHO] * sqrt(x1); // coefficient of shear viscosity = (kinematic viscosity)*rho
  *nu2 = 0.0; // coefficient of bulk viscosity
}
