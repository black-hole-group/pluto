Torus 2D hydro simulation, Newtonian
========================================

Two-dimensional (axially symmetric) hydrodynamical simulation of a RIAF without cooling. The calculations begin from an equilibrium configuration consisting of a thick torus (Papaloizou-Pringle) with constant specific angular momentum. Accretion is induced by the addition of a small anomalous azimuthal shear stress which. 

For a tutorial on how to run this simulation, please visit [this page](https://github.com/black-hole-group/group-wiki/wiki/Pluto).

# Code units

$$G=M=1$$

We usually parametrize time in units of the orbital time at the radius of maximum density (or pressure), 
$$t_{\rm orb} (r_0) = \frac{2 \pi r_0}{v_\phi(r0)} $$
where $v_\phi(r0)$ is the Keplerian velocity for the appropriate potential at $r_0$. 

# Viscous stress tensor

Parametrized following Stone et al. (1999) with the kinematic viscosity given by 
$$\nu = \alpha \frac{c_s^2}{\Omega_k} \approx \alpha \rho$$
The shear stress is defined in file `visc_nu.c`.

# Gravity

A Newtonian potential is adopted, 
$$\Phi = - \frac{GM}{r}$$ 

The potential is defined in `init.c`.

# Coordinate system and grid

Coordinate system: spherical polar

Resolution: two different resolutions are suggested, one for laptop tests in `pluto_laptop.ini` and another for cluster tests in `pluto_cluster.ini`. The radial range changes accordingly.

## Non-uniform grid in $\theta$

If desired, a non-uniform grid in the $\theta$-direction can be adopted and is specified in `pluto_theta.ini`. 

The reasoning is that we are not necessarily interested in having resolution towards the poles, since we are not simulating relativistic jets. According to Moscibrodzka et al. (2014), the jet is produced inside a cone of angle $25\circ$. Therefore, compared to a uniform grid in $\theta$, our approach is to have 1/3 of the uniform grid resolution within the "jet cone".

The formula for calculating the required number of elements inside one half of the jet cone is the following:
$$N_\theta ({\rm jet}) = N_{\theta} ({\rm total}) \left( \frac{25^\circ}{180^\circ} \right) ({\rm reduction \ factor})$$
where we typically use `reduction factor = 1/3`.

For a concrete example, suppose we choose $N_\theta=100$ and decide to decrease the resolution to 1/3 within $25\circ$ of the poles. According to the formula above, one has $n=5$ inside each half of the cone. Therefore, we would have 5 elements in the upper half of the cone, 5 in the second half and 90 elements in the rest of the theta domain. 

In Pluto, you can insert this modified resolution in `pluto.ini` as follows:

    X2-grid  3    0.0  5 u 0.436 90 u 2.705 5 u 3.1415



