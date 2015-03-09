Integrators 
============


Integrators in libmd               {#md-integrator}
=====================

In libmd there are integrators for flat space (discussed here) and integrators for curved space. 
All integrators in libmd are symplectic, see [1](#md-int-ref) and [2](#md-int-ref). 
Symplectic integrators are integrators that preserve the Lagrangian symmetries.
Therefore integration methods like explicit Euler or Runge-Kutta are not used.

Flat space integrators                     {#md-integratorflat}
-------------------

For flat space molecular dynamics two integrators are provided: 

- [Symplectic Euler](#md-symeul) (also called semi-implicit Euler, semi-explicit Euler, Euler–Cromer, and Newton–Størmer–Verlet)
- [Velocity Verlet](#md-vverlet) 

These 

### Symplectic Euler       {#md-symeul}
Symplectic Euler is a first order method of the form
\f{align*}{
    \dot{x}^{\mu}_{i+1} &= \dot{x}^{\mu}_{i}+\tfrac{h}{m_i}F^{\mu}_{i} \\
    x^{\mu}_{i+1} &= x^{\mu}_{i}+h\dot{x}^{\mu}_{i+1}
\f}
where \f$h\f$ is a timestep. 

### Velocity Verlet        {#md-vverlet}
Velocity Verlet is a second order method of the form
\f{align}{
    x^{\mu}_{i+1} &= x^{\mu}_{i} + h \dot{x}^{\mu}_{i} + \tfrac{h^2}{2 m} F^{\mu}_{i} \\
    \dot{x}^{\mu}_{i+1} &= \dot{x}^{\mu}_{i} + \tfrac{h^2}{2 m} (F^{\mu}_{i}+F^{\mu}_{i+1})
\f}
where \f$h\f$ is a time step. 

Integrators structure                     {#md-integrators}
-------------------
Integrators in libmd can be controlled with the integrate structure in <tt>md<dim>::integrator</tt> of type <tt>#integrators</tt>.
For flat space two objects are relevant:
* <tt>integrators::h</tt>      the time step
* <tt>integrators::method</tt>     the integration method

The time step controls the time evolution.
As a rule of thumb the time steps shouldn't be to big because the velocities will arbitrarily rise as the potentials give too large energies, as particles end up where they shouldn't be (due to errors). 
A time step should not be to large as floating point errors will accumulate creating an error and/or barely any evolution takes place.
Keep in mind that the physics should not be dependent on the choice of time step, it it does, than the time step is probably too small or too large. 
If the time step is not to small and not to large the symplecticity of the integrator will ensure the Lagrangian symmetries are preserved. 
The default time step size is \f$h=10^{-3}\f$.

The flat space integration methods supported by libmd are enumerated in #INTEGRATOR::integrator.
The default method is velocity Verlet.

### References             {#md-int-ref}
1. [Discrete mechanics and variational integrators JE Marsden, M West - Acta Numerica 2001, 2001 - Cambridge Univ Press](http://dx.doi.org/10.1017/S096249290100006X)
2. [Discrete geometric mechanics for variational time integrator nt-refs A Stern, M Desbrun - ACM SIGGRAPH 2006 Courses, 2006](http://dx.doi.org/10.1145/1185657.1185669)

