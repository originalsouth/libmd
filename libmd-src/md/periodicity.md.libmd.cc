#define __libmd_src_file__
#ifndef libmd_h
#include "../../libmd.h"
#endif

template<ui dim> void md<dim>::thread_periodicity(ui i)
{
    //!
    //! This function loops through the <tt>dim</tt> dimensions, and updates
    //! the position and velocity of particle <tt>i</tt> to respect the boundary
    //! conditions of any boundary it might have passed through in the last time step.
    //!
    for(ui d=0;d<dim;d++) boundary(d,i,this);
}

template<ui dim> void md<dim>::thread_periodicity(ldf x[dim])
{
    //!
    //! This function loops through the <tt>dim</tt> dimensions, and updates
    //! the position and velocity of particle <tt>i</tt> to respect the boundary
    //! conditions of any boundary it might have passed through in the last time step.
    //!
    for(ui d=0;d<dim;d++) boundary(d,x,this);
}

template<ui dim> void md<dim>::periodicity()
{
    //!
    //! Update positions and velocities of particles to respect the appropriate boundary conditions.
    //! <br>
    //! This function loops through the <tt>dim</tt> dimensions, and updates the positions
    //! and velocities of all (non-fixed) particles that have passed through the boundary perpendicular
    //! to dimension <tt>d</tt>,
    //! based on the value of <tt>simbox.bcond[d]</tt>.
    //!
    for(ui i=0;i<N;i++) if(!particles[i].fix)
    {
        if(particles[i].usepbcond) for(ui d=0;d<dim;d++) boundary(particles[i].pbc[d],d,i,this);
        else for(ui d=0;d<dim;d++) boundary(d,i,this);
    }
}

template<ui dim> void md<dim>::update_boundaries()
{
    //!
    //! Update the box matrix using the shear velocities.
    //! Increments the box matrix \f$L_{ij}\f$ (stored in <tt>simbox.Lshear</tt>)
    //! by \f$h\times v_{ij}\f$ where \f$v_{ij}\f$ is the shear velocity matrix
    //! (stored in <tt>simbox.vshear</tt>) and \f$h\f$ is the integration timestep
    //! (stored in <tt>integrator.h</tt>).
    //! <br>
    //! Also shifts \f$L_{ij}\f$ to ensure that it remains within the bounds \f$-L_{ii}/2 \leq L_{ij} \leq L_{ii}/2\f$.
    //!
    // update box matrix for shear
    for(ui j=0;j<dim;j++) for (ui k=0; k<dim; k++)
    {
        simbox.Lshear[j][k] += simbox.vshear[j][k]*integrator.h;
        // shift by appropriate box lengths so that the off-diagonal entries are in the range -L[j][j]/2 to L[j][j]/2 consistent with the positions
        if (j != k && (simbox.bcond[j]==BCOND::PERIODIC || simbox.bcond[j]==BCOND::BOXSHEAR))
        {
            while(simbox.Lshear[j][k]>simbox.Lshear[j][j]/2.) simbox.Lshear[j][k]-=simbox.Lshear[j][j];
            while(simbox.Lshear[j][k]<-simbox.Lshear[j][j]/2.) simbox.Lshear[j][k]+=simbox.Lshear[j][j];
        }
    }
    simbox.invert_box();
}
