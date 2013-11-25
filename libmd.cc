///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  Begin of LIBRARY source                                                                                      //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#include "libmd.h"

#include "libmd/potentials.libmd.cc"
#include "libmd/particle.libmd.cc"
#include "libmd/interact.libmd.cc"
#include "libmd/pairpotentials.libmd.cc"
#include "libmd/integrators.libmd.cc"

template<ui dim> md<dim>::md()
{
    N=0;
}

template<ui dim> md<dim>::md(ui particlenr)
{
    N=particlenr;
    particles.resize(N);
    network.skins.resize(N);
}

template<ui dim> void md<dim>::add_typeinteraction(ui type1,ui type2,ui potential,vector<ldf> *parameters)
{
    pair<ui,ui> id=network.hash(type1,type2);
    interactiontype itype(potential,parameters,v(potential,network.rco,network.rcosq,parameters));
    if(network.lookup.find(id)==network.lookup.end()) network.library.push_back(itype),network.lookup[id]=network.library.size()-1;
    else network.library[network.lookup[id]]=itype;
}

template<ui dim> void md<dim>::mod_typeinteraction(ui type1,ui type2,ui potential,vector<ldf> *parameters)
{
    pair<ui,ui> id=network.hash(type1,type2);
    interactiontype itype(potential,parameters,v(potential,network.rco,network.rcosq,parameters));
    if(network.lookup.find(id)==network.lookup.end()) network.library.push_back(itype),network.lookup[id]=network.library.size()-1;
    else network.library[network.lookup[id]]=itype;
}

template<ui dim> ldf md<dim>::distsq(ui p1,ui p2)
{
    ldf retval=0.0;
    for(ui i=0;i<dim;i++)
    {
        ldf ad=fabs(particles[p2].x[i]-particles[p1].x[i]),d;
        switch(simbox.bcond[i])
        {
            case 1: d=(ad<simbox.L[i]/2.0?ad:simbox.L[i]-ad); break;
            default: d=ad; break;
        }
        retval+=pow(d,2);
    }
    return retval;
}

template<ui dim> ldf md<dim>::dd(ui i,ui p1,ui p2)
{
    ldf ad=particles[p2].x[i]-particles[p1].x[i],d;
    switch(simbox.bcond[i])
    {
        case 1: d=fmod(3.0*simbox.L[i]/2.0+ad,simbox.L[i])-simbox.L[i]/2.0; break;
        default: d=ad; break;
    }
    return d;
}

template<ui dim> void md<dim>::thread_clear_forces(ui i)
{
    for(ui d=0;d<dim;d++) particles[i].F[d]=0.0;
}

//TODO: Implement atomics
//TODO: What if potential is velocity dependent (damping)?
template<ui dim> void md<dim>::thread_calc_forces(ui i)
{
    for(ui j=network.skins[i].size()-1;j<numeric_limits<ui>::max();j--) if(i>network.skins[i][j].neighbor)
    {
        const ldf rsq=distsq(i,network.skins[i][j].neighbor);
        if(rsq<network.rcosq)
        {
            const ldf r=sqrt(rsq);
            const ldf dVdr=v.dr(network.library[network.skins[i][j].interaction].potential,r,rsq,&network.library[network.skins[i][j].interaction].parameters);
            for(ui d=0;d<dim;d++)
            {
                const ldf delta=dd(d,i,network.skins[i][j].neighbor);
                const ldf F=delta*dVdr/r;
                particles[i].F[d]+=F; //FIXME: atomic
                particles[network.skins[i][j].neighbor].F[d]-=F; //FIXME: atomic
            }
        }
    }
}

//FIXME: Thomas
template<ui dim> void md<dim>::index()
{
    for(ui i=0;i<N;i++)
    {
        network.skins[i].clear();
        for(ui j=0;j<N;j++) if(i!=j and distsq(i,j)<network.sszsq)
        {
            interactionneighbor in(j,network.lookup[network.hash(i,j)]);
            network.skins[i].push_back(in);
        }
    }
}

//TODO: Make parallel launcher
template<ui dim> void md<dim>::calc_forces()
{
    for(ui i=0;i<N;i++) thread_clear_forces(i);
    for(ui i=0;i<N;i++) thread_calc_forces(i);
}

//TODO: Make parallel launcher
template<ui dim> void md<dim>::recalc_forces()
{
    for(ui i=0;i<N;i++) thread_calc_forces(i);
}

template<ui dim> void md<dim>::thread_integrate(ui i,ui gen)
{
    switch(integrator.method)
    {
        case 1:
            switch(gen)
            {
                case 0:
                    for(ui d=0;d<dim;d++)
                    {
                        particles[i].xp[d]=particles[i].x[d];
                        particles[i].x[d]+=integrator.h*particles[i].dx[d]+0.5*pow(integrator.h,2)*particles[i].F[d];
                    }
                break;
                case 1:
                    for(ui d=0;d<dim;d++) particles[i].dx[d]+=0.5*integrator.h*particles[i].F[d];
                break;
            }
        break;
        default:
            const ldf o=integrator.h/particles[i].m;
            for(ui d=0;d<dim;d++)
            {
                particles[i].dx[d]+=o*particles[i].F[d];
                particles[i].xp[d]=particles[i].x[d];
                particles[i].x[d]+=integrator.h*particles[i].dx[d];
            }
        break;
    }

}

//TODO: Make parallel launcher
//TODO: Implement boundary conditions
//TODO: Implement masses
template<ui dim> void md<dim>::integrate()
{
    switch(integrator.method)
    {
        case 1:
        for(ui i=0;i<N;i++) thread_integrate(i,0);
        recalc_forces();
        for(ui i=0;i<N;i++) thread_integrate(i,1);
        break;
        default:
        for(ui i=0;i<N;i++) thread_integrate(i,0);
        break;
    }


}

template<ui dim> void md<dim>::timestep()
{
    if(network.update) index();
    calc_forces();
    integrate();
}

template<ui dim> void md<dim>::timesteps(ui k)
{
    for(ui i=0;i<k;i++) timestep();
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  End of LIBRARY source                                                                                        //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
