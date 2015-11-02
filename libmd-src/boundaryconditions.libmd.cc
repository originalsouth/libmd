#define __libmd_src_file__
#ifndef libmd_h
#include "../libmd.h"
#endif

template<ui dim> void BCOND_NONE(ui d,ui i,void *sys)
{
    //!
    //! Periodicity function to be called if
    //! dimension <tt>d</tt> has no boundary conditions.
    //!
    //! This function does nothing
    //!
    (void) d;
    (void) i;
    (void) sys;
}

template<ui dim> void BCOND_NONE(ui d,ldf x[dim],void *sys)
{
    //!
    //! Periodicity function to be called if
    //! dimension <tt>d</tt> has no boundary conditions.
    //!
    //! This function does nothing
    //!
    (void) d;
    (void) x;
    (void) sys;
}

template<ui dim> void BCOND_PERIODIC(ui d,ui i,void *sys)
{
    //!
    //! Periodicity function to be called if
    //! dimension <tt>d</tt> has periodic boundary conditions.
    //!
    //! Checks if particle <tt>i</tt> has crossed the boundary perpendicular to
    //! dimension <tt>d</tt> and, if so, shifts its coordinate in that dimension by multiples of
    //! <tt>simbox.L[d]</tt> so that it is within the bounds <tt>(-simbox.L[d]/2,simbox.L[d]/2)</tt>.
    //!
    using namespace std;
    const ldf dx=SYS->simbox.L[d]*round(SYS->particles[i].x[d]/SYS->simbox.L[d]);
    SYS->particles[i].xp[d]-=dx;
    SYS->particles[i].x[d]-=dx;
}

template<ui dim> void BCOND_PERIODIC(ui d,ldf x[dim],void *sys)
{
    //!
    //! Periodicity function to be called if
    //! dimension <tt>d</tt> has periodic boundary conditions.
    //!
    //! Checks if point <tt>x</tt> is outside the boundary perpendicular to
    //! dimension <tt>d</tt> and, if so, shifts its coordinate in that dimension by multiples of
    //! <tt>simbox.L[d]</tt> so that it is within the bounds <tt>(-simbox.L[d]/2,simbox.L[d]/2)</tt>.
    //!
    using namespace std;
    const ldf dx=SYS->simbox.L[d]*round(x[d]/SYS->simbox.L[d]);
    x[d]-=dx;
}

template<ui dim> void BCOND_HARD(ui d,ui i,void *sys)
{
    //!
    //! Periodicity function to be called if
    //! dimension <tt>d</tt> has hard boundary conditions.
    //!
    //! Checks if particle <tt>i</tt> has crossed the boundary perpendicular to
    //! dimension <tt>d</tt> and, if so, updates its position and velocity
    //! to respect a hard wall reflection. The particle position is mirrored
    //! across the boundary wall, whereas its velocity component perpendicular
    //! to the boundary wall is reversed.
    //! <br>
    //! This function correctly takes into account skewed boundary conditions,
    //! and uses the box matrices <tt>simbox.Lshear</tt> and <tt>simbox.vshear</tt>
    //! to calculate the reflections if <tt>simbox.useLshear</tt> is <tt>true</tt>.
    //!
    using namespace std;
    if (SYS->simbox.useLshear)
    {
        ldf s=0;
        for (ui k=0;k<dim;k++) s+=SYS->simbox.LshearInv[d][k]*SYS->particles[i].x[k];
        if (abs(s) > 0.5) // particle has hit the hard boundary as distorted by the shear
        {
            if (abs(s) > 1.) { WARNING("dynamics led to particle displacement bigger than box size; hard boundary reflections undefined"); }
            ldf nhat[dim];
            ldf nlen=0.,vperp=0.,xperp=0.,x0perp;

            // the normal vector to the box boundary in dimension d is the dth row of LshearInv
            for (ui k=0;k<dim;k++) nlen += SYS->simbox.LshearInv[d][k]*SYS->simbox.LshearInv[d][k];
            nlen = sqrt(nlen);

            // projection of velocity and position perpendicular to boundary wall
            for (ui k=0;k<dim;k++) {
                nhat[k] = SYS->simbox.LshearInv[d][k]/nlen;
                vperp += nhat[k]*SYS->particles[i].dx[k];
                xperp += nhat[k]*SYS->particles[i].x[k];
            }

            x0perp = nhat[d]*SYS->simbox.Lshear[d][d]*0.5*(s > 0.? 1.:-1.);

            // subtract perpendicular component twice to get reflected velocity
            for (ui k=0;k<dim;k++)
            {
                SYS->particles[i].dx[k] -= 2.0*vperp*nhat[k];
                SYS->particles[i].x[k] -= 2.0*(xperp-x0perp)*nhat[k]; // reflection about a plane passing through point set by x0perp
                SYS->particles[i].xp[k] += 2.0*(xperp-x0perp)*nhat[k];
            }
        }
    }
    else
    {
        const ldf xnew=SYS->simbox.L[d]*(abs(SYS->particles[i].x[d]/SYS->simbox.L[d]+0.5-2.0*floor(SYS->particles[i].x[d]/(2.0*SYS->simbox.L[d])+0.75))-0.5);
        const ldf sign=(((int)round(SYS->particles[i].x[d]/SYS->simbox.L[d]))%2?-1.0:1.0);
        SYS->particles[i].xp[d]+=sign*(xnew-SYS->particles[i].x[d]);
        SYS->particles[i].x[d]=xnew;
        SYS->particles[i].dx[d]*=sign;
    }
}

template<ui dim> void BCOND_HARD(ui d,ldf x[dim],void *sys)
{
    //!
    //! Periodicity function to be called if
    //! dimension <tt>d</tt> has hard boundary conditions.
    //!
    //! Checks if point <tt>x</tt> is outside the boundary perpendicular to
    //! dimension <tt>d</tt> and, if so, updates its position
    //! to respect a hard wall reflection. The particle position is mirrored
    //! across the boundary wall.
    //! <br>
    //! This function correctly takes into account skewed boundary conditions,
    //! and uses the box matrices <tt>simbox.Lshear</tt> and <tt>simbox.vshear</tt>
    //! to calculate the reflections if <tt>simbox.useLshear</tt> is <tt>true</tt>.
    //!
    using namespace std;
    if (SYS->simbox.useLshear)
    {
        ldf s=0;
        for (ui k=0;k<dim;k++) s+=SYS->simbox.LshearInv[d][k]*x[k];
        if (abs(s) > 0.5) // particle has hit the hard boundary as distorted by the shear
        {
            if (abs(s) > 1.) { WARNING("dynamics led to particle displacement bigger than box size; hard boundary reflections undefined"); }
            ldf nhat[dim];
            ldf nlen=0.,xperp=0.,x0perp;
            // the normal vector to the box boundary in dimension d is the dth row of LshearInv
            for (ui k=0;k<dim;k++) nlen += SYS->simbox.LshearInv[d][k]*SYS->simbox.LshearInv[d][k];
            nlen = sqrt(nlen);
            // projection of velocity and position perpendicular to boundary wall
            for (ui k=0;k<dim;k++) {
                nhat[k] = SYS->simbox.LshearInv[d][k]/nlen;
                xperp += nhat[k]*x[k];
            }
            x0perp = nhat[d]*SYS->simbox.Lshear[d][d]*0.5*(s > 0.? 1.:-1.);

            // subtract perpendicular component twice to get reflected velocity
            for (ui k=0;k<dim;k++) x[k] -= 2.0*(xperp-x0perp)*nhat[k]; // reflection about a plane passing through point set by x0perp
        }
    }
    else
    {
        const ldf xnew=SYS->simbox.L[d]*(abs(x[d]/SYS->simbox.L[d]+0.5-2.0*floor(x[d]/(2.0*SYS->simbox.L[d])+0.75))-0.5);
        x[d]=xnew;
    }
}

template<ui dim> void BCOND_BOXSHEAR(ui d,ui i,void *sys)
{
    //!
    //! Periodicity function to be called if
    //! dimension <tt>d</tt> has sheared boundary conditions.
    //!
    //! Checks if particle <tt>i</tt> has crossed the boundary perpendicular to
    //! dimension <tt>d</tt> and, if so, updates its position and velocity
    //! according to the box shear matrices stored in <tt>simbox.Lshear</tt>
    //! and <tt>simbox.vshear</tt>. The particle position
    //!
    using namespace std;
    const ldf boundaryCrossing=round(SYS->particles[i].x[d]/SYS->simbox.L[d]);
    if((int)boundaryCrossing) for(ui k=0;k<dim;k++)
    {
        SYS->particles[i].x[k]-=SYS->simbox.Lshear[k][d]*boundaryCrossing;
        SYS->particles[i].xp[k]-=SYS->simbox.Lshear[k][d]*boundaryCrossing;
        SYS->particles[i].dx[k]-=SYS->simbox.vshear[k][d]*boundaryCrossing;
    }
}

template<ui dim> void BCOND_BOXSHEAR(ui d,ldf x[dim],void *sys)
{
    //!
    //! Periodicity function to be called if
    //! dimension <tt>d</tt> has sheared boundary conditions.
    //!
    //! Checks if point <tt>x</tt> is outside the boundary perpendicular to
    //! dimension <tt>d</tt> and, if so, updates its position
    //! according to the box shear matrices stored in <tt>simbox.Lshear</tt>
    //! and <tt>simbox.vshear</tt>. The particle position
    //!
    using namespace std;
    const ldf boundaryCrossing=round(x[d]/SYS->simbox.L[d]);
    if((int)boundaryCrossing) for(ui k=0;k<dim;k++) x[k]-=SYS->simbox.Lshear[k][d]*boundaryCrossing;
}
