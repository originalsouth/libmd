#ifndef libmd_h
#include "../../libmd.h"
#endif


template<ui dim> void md<dim>::thread_periodicity_periodic(ui d,ui i)
{
    ldf dx=simbox.L[d]*round(particles[i].x[d]/simbox.L[d]);
    particles[i].xp[d]-=dx;
    particles[i].x[d]-=dx;
}


template<ui dim> void md<dim>::thread_periodicity_boxshear(ui d,ui i) 
{
    ldf boundaryCrossing=round(particles[i].x[d]/simbox.L[d]);
    if(fabs(boundaryCrossing)>0.1) for(ui k=0;k<dim;k++)
    {
        particles[i].x[k]-=simbox.Lshear[k][d]*boundaryCrossing;
        particles[i].xp[k]-=simbox.Lshear[k][d]*boundaryCrossing;
        particles[i].dx[k]-=simbox.vshear[k][d]*boundaryCrossing;
    }
}

template<ui dim> void md<dim>::thread_periodicity_hard(ui d,ui i)
{   
    if (simbox.boxShear)
    {
        ldf s=0;
        for (ui k=0;k<dim;k++) s+=simbox.LshearInv[d][k]*particles[i].x[k]; 
        if (fabs(s) > 0.5) // particle has hit the hard boundary as distorted by the shear
        {
            if (fabs(s) > 1.) { WARNING("Dynamics led to particle displacement bigger than box size. Hard boundary reflections undefined."); }  
            ldf nhat[dim];
            ldf nlen=0.,vperp=0.,xperp=0.,x0perp;
            
            // the normal vector to the box boundary in dimension d is the dth row of LshearInv
            for (ui k=0;k<dim;k++) nlen += simbox.LshearInv[d][k]*simbox.LshearInv[d][k];
            nlen = sqrt(nlen);
            
            // projection of velocity and position perpendicular to boundary wall
            for (ui k=0;k<dim;k++) {
                nhat[k] = simbox.LshearInv[d][k]/nlen;
                vperp += nhat[k]*particles[i].dx[k]; 
                xperp += nhat[k]*particles[i].x[k]; 
            }
            
            x0perp = nhat[d]*simbox.Lshear[d][d]*0.5*(s > 0.? 1.:-1.);
            
            // subtract perpendicular component twice to get reflected velocity
            for (ui k=0;k<dim;k++) 
            { 
                particles[i].dx[k] -= 2*vperp*nhat[k];
                particles[i].x[k] -= 2*(xperp-x0perp)*nhat[k]; // reflection about a plane passing through point set by x0perp
                particles[i].xp[k] += 2*(xperp-x0perp)*nhat[k]; 
            }
        }
    }
    else 
    {
        ldf xnew=simbox.L[d]*(fabs(particles[i].x[d]/simbox.L[d]+0.5-2.0*floor(particles[i].x[d]/(2.0*simbox.L[d])+0.75))-0.5);
        ldf sign=(((int)round(particles[i].x[d]/simbox.L[d]))&1?-1.0:1.0);
        particles[i].xp[d]+=sign*(xnew-particles[i].x[d]);
        particles[i].x[d]=xnew;
        particles[i].dx[d]*=sign;
    }
}

template<ui dim> void md<dim>::thread_periodicity(ui i)
{
    if(simbox.bcond) for(ui d=0;d<dim;d++) switch(simbox.bcond[d])
    {
        case BCOND::PERIODIC: if(!particles[i].fix) thread_periodicity_periodic(d,i); break;
        case BCOND::BOXSHEAR: if(!particles[i].fix) thread_periodicity_boxshear(d,i); break;
        case BCOND::HARD: if(!particles[i].fix) thread_periodicity_hard(d,i); break;
    }
}

template<ui dim> void md<dim>::periodicity()
{
    if(simbox.bcond) for(ui d=0;d<dim;d++) switch(simbox.bcond[d])
    {
        case BCOND::PERIODIC:
        {
            DEBUG_2("applying periodic boundary conditions in %u dimension",d);
            for(ui i=0;i<N;i++) if(!particles[i].fix) thread_periodicity_periodic(d,i);
        }
        break;
        case BCOND::BOXSHEAR:
        {
            DEBUG_2("applying boxshear boundary conditions in %u dimension",d);
            for(ui i=0;i<N;i++) if(!particles[i].fix) thread_periodicity_boxshear(d,i);
        }
        break;
        case BCOND::HARD:
        {
            DEBUG_2("applying hard boundary conditions in %u dimension",d);
            for(ui i=0;i<N;i++) if(!particles[i].fix) thread_periodicity_hard(d,i);
        }
        break;
    }
}


template<ui dim> void md<dim>::update_boundaries()
{
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

