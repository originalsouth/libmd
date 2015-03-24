Monge patch integrators 
============


Integrators in libmd               {#mp-integrator}
=====================

In libmd there are integrators for [flat space](#md-integrator) and integrators for curved space (discussed here). 
All integrators in libmd are symplectic, see [1](#md-int-ref) and [2](#md-int-ref). 
Symplectic integrators are integrators that preserve the Lagrangian symmetries.
Therefore integration methods like explicit Euler or Runge-Kutta are not used.

Monge patch space integrators                     {#mp-integrator-patch}
-------------------

For curved space molecular dynamics multiple integrators a provide some of which are invalid:
- [Symplectic Euler](@ref md-symeul) inherited from flat space
- [Velocity Verlet](@ref md-vverlet) inherited from flat space
- [van Zuiden](#mp-vanzuiden) possibly somewhat narcisstically named.
- [van Zuiden WFI](#mp-vanzuiden) the van Zuiden integrator without fixed point iterations
- [van Zuiden P](#mp-vanzuiden) the van Zuiden integrator with limited fixed point iterations

### The van Zuiden integrator       {#mp-vanzuiden}
The van Zuiden integrator is derived using Variational integrators -- see [1](@ref md-int-ref) and [2](@ref md-int-ref) -- using the following Lagrangian:
\f{align}{
    L=\tfrac{1}{2} m g_{\mu \nu} \dot{x}^{\mu} \dot{x}^{\nu} - V(x^{\rho})
\f}
where \f$g_{\mu \nu}\f$ is the metric tensor.
By discritizing in time and applying the discritized Euler--Lagrange equatioins we obtain:
\f{align}{
    \epsilon^{\rho}=g^{\sigma \rho} \left( \Gamma_{\nu \sigma \mu} \epsilon^{\mu} \epsilon^{\nu} + C_{\sigma} \right)
\f}
where \f$g^{\mu \nu}\f$ is the metric inverse, \f$\epsilon^{\mu}=x^{\mu}_{t+1}-x^{\mu}_{t}\f$, \f$\Gamma_{\sigma \mu \nu}\f$ are the Christoffel symbols of the first kind.
The latter equations can be solved using fixed points iterations using the starting point:
\f{align}{
    \epsilon^{\rho}=g^{\sigma \rho} C_{\sigma}
\f}

The metric can be modified as discussed [here](@ref mp-curvature).

Integrators structure                     {#mp-integrators}
-------------------
Integrators in libmd can be controlled with the integrate structure in <tt>md<dim>::integrator</tt> of type <tt>#integrators</tt>.
For curved space two objects are relevant:
* <tt>integrators::h</tt>      the time step
* <tt>integrators::method</tt>     the integration method
* <tt>integrators::generations</tt>     number of allowed fixed point iterations

The time step controls the time evolution.
As a rule of thumb the time steps shouldn't be to big because the velocities will arbitrarily rise as the potentials give too large energies, as particles end up where they shouldn't be (due to errors). 
A time step should not be to large as floating point errors will accumulate creating an error and/or barely any evolution takes place.
Keep in mind that the physics should not be dependent on the choice of time step, it it does, than the time step is probably too small or too large. 
If the time step is not to small and not to large the symplecticity of the integrator will ensure the Lagrangian symmetries are preserved. 
The default time step size is \f$h=10^{-3}\f$. The default value of generations is 10.

The curved space integration methods supported by libmd are enumerated in #MP_INTEGRATOR::mp_integrator.
The default method is van Zuiden.
