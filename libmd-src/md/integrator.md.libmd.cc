#define __libmd_src_file__
#ifndef libmd_h
#include "../../libmd.h"
#endif

template<ui dim> void md<dim>::thread_seuler(ui i)
{
    //!
    //! This function integrates particle <tt>i</tt> using symplectic Euler. <br>
    //! The symplectic Euler integrator is simpelest (and fastest) symplectic integrator. <br>
    //! The integrator is a first order integrator. <br>
    //! Sometimes symplectic Euler is also refered: to as Semi-implicit Euler, Euler--Cromer or Newton-Stoermer-Verlet. <br>
    //! Symplectic Euler takes the following form:
    //! \f{align}{ \dot{x}^{\mu}_{t+1}=&\dot{x}^{\mu}_{t}+\tfrac{h}{m} F^{\mu}_{t}\\ x^{\mu}_{t+1}=& x^{\mu}_{t}+h\dot{x}^{\mu}_{t+1}\f}
    //!
    ldf o=integrator.h/particles[i].m;
    memcpy(particles[i].xp,particles[i].x,dim*sizeof(ldf));
    for(ui d=0;d<dim;d++)
    {
        particles[i].dx[d]+=o*particles[i].F[d];
        particles[i].x[d]+=integrator.h*particles[i].dx[d];
    }
}

template<ui dim> void md<dim>::thread_vverlet_x(ui i)
{
    //!
    //! This function integrates particle position <tt>i</tt> using Velocity Verlet. <br>
    //! See md<dim>::thread_vverlet_dx for the particle velocity part. <br>
    //! Velocity Verlet is a second order symplectic integrator. <br>
    //! It is slower than symplectic Euler (see md<dim>::thread_seuler), yet more accurate. <br>
    //! Velocity Verlet takes the following form:
    //! \f{align}{x^{\mu}_{t+1}=&x^{\mu}_{t}+h\dot{x}^{\mu}_{t}+\tfrac{h^2}{2}F^{\mu}_{t}\\ \dot{x}^{\mu}_{t+1}=&\dot{x}^{\mu}_{t}+\tfrac{h}{2}(F^{\mu}_{t}+F^{\mu}_{t+1})\f}
    //! Because of the second force (\f$F^{\mu}_{t+1}\f$) calculation this function is split in two. <br>
    //!
    memcpy(particles[i].xp,particles[i].x,dim*sizeof(ldf));
    for(ui d=0;d<dim;d++) particles[i].x[d]+=integrator.h*particles[i].dx[d]+0.5*std::pow(integrator.h,2)*particles[i].F[d]/particles[i].m;
}

template<ui dim> void md<dim>::thread_vverlet_dx(ui i)
{
    //!
    //! This function integrates particle velocity <tt>i</tt> using Velocity Verlet. <br>
    //! See md<dim>::thread_vverlet_x for the particle position part. <br>
    //! Velocity Verlet is a second order symplectic integrator. <br>
    //! It is slower than symplectic Euler (see md<dim>::thread_seuler), yet more accurate. <br>
    //! Velocity Verlet takes the following form:
    //! \f{align}{x^{\mu}_{t+1}=&x^{\mu}_{t}+h\dot{x}^{\mu}_{t}+\tfrac{h^2}{2}F^{\mu}_{t}\\ \dot{x}^{\mu}_{t+1}=&\dot{x}^{\mu}_{t}+\tfrac{h}{2}(F^{\mu}_{t}+F^{\mu}_{t+1})\f}
    //! Because of the second force (\f$F^{\mu}_{t+1}\f$) calculation this function is split in two. <br>
    //!
    for(ui d=0;d<dim;d++) particles[i].dx[d]+=0.5*integrator.h*particles[i].F[d]/particles[i].m;
}

template<ui dim> void md<dim>::thread_first_order(ui i)
{
    //!
    //! This function integrates particle position <tt>i</tt> using Euler first order. <br>
    //! This is useful for implementing Vicsek type models
    //!
    memcpy(particles[i].xp,particles[i].x,dim*sizeof(ldf));
    for(ui d=0;d<dim;d++) particles[i].x[d]+=integrator.h*particles[i].dx[d];
}

template<ui dim> void md<dim>::thread_overdamped(ui i)
{
    //!
    //! This function integrates particle position <tt>i</tt> using Euler first order. <br>
    //! This is useful for implementing overdamped systems
    //!
    memcpy(particles[i].xp,particles[i].x,dim*sizeof(ldf));
    for(ui d=0;d<dim;d++) particles[i].dx[d]=-particles[i].F[d]/avars.overdamped_gamma;
    for(ui d=0;d<dim;d++) particles[i].x[d]+=integrator.h*particles[i].dx[d];
}

template<ui dim> void md<dim>::integrate()
{
    //!
    //! This function calls one of the specified available integrators depending on the value of md<dim>::integrator.method. <br>
    //! The default call is symplectic Euler (see md<dim>::thread_seuler).
    //!
    avars.export_force_calc=true;
    switch(integrator.method)
    {
        case INTEGRATOR::FO_OVERDAMPED:
            DEBUG_2("integrating using overdamped first order (Euler)");
            for(ui i=0;i<N;i++) if(!particles[i].fix) thread_overdamped(i);
        break;
        case INTEGRATOR::FO:
            DEBUG_2("integrating using first order (Euler)");
            for(ui i=0;i<N;i++) if(!particles[i].fix) thread_first_order(i);
        break;
        case INTEGRATOR::VVERLET:
            DEBUG_2("integrating using symplectic Velocity Verlet");
            for(ui i=0;i<N;i++) if(!particles[i].fix) thread_vverlet_x(i);
            recalc_forces();
            for(ui i=0;i<N;i++) if(!particles[i].fix) thread_vverlet_dx(i);
        break;
        default:
            DEBUG_2("integrating using symplectic Euler method");
            for(ui i=0;i<N;i++) if(!particles[i].fix) thread_seuler(i);
        break;
    }
    periodicity();
}

template<ui dim> void md<dim>::timestep()
{
    //!
    //! This function defines a timestep by combining the md<dim>::calc_forces and md<dim>::integrate functions. <br>
    //! If there is a boxshear present it updates the boundaries.
    //!
    if(simbox.useLshear)
    {
        DEBUG_2("updating boxshear boundaries");
        update_boundaries();
    }
    calc_forces();
    integrate();
    run_hooks();
}

template<ui dim> void md<dim>::timesteps(ui k)
{
    //!
    //! This function calls md<dim>::timestep <tt>k</tt> times.<br>
    //! (The looped is sandwiched between DEBUG_TIMER statements.)
    //!
    DEBUG_TIMER("start timesteps");
    for(ui i=0;i<k;i++) timestep();
    DEBUG_TIMER("stop timesteps");
}
