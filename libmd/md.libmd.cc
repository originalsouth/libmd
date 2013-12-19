#ifndef libmd_h
#include "../libmd.h"
#endif

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

template<ui dim> bool md<dim>::add_typeinteraction(ui type1,ui type2,ui potential,vector<ldf> *parameters)
{
    pair<ui,ui> id=network.hash(type1,type2);
    interactiontype itype(potential,parameters,v(potential,network.rco,parameters));
    if(network.lookup.find(id)==network.lookup.end())
    {
        network.library.push_back(itype),network.lookup[id]=network.library.size()-1;
        network.backdoor.push_back(id);
        return true;
    }
    else return false;
}

template<ui dim> bool md<dim>::mod_typeinteraction(ui type1,ui type2,ui potential,vector<ldf> *parameters)
{
    pair<ui,ui> id=network.hash(type1,type2);
    interactiontype itype(potential,parameters,v(potential,network.rco,network.rcosq,parameters));
    if(network.lookup.find(id)==network.lookup.end()) return false;
    else
    {
        network.library[network.lookup[id]]=itype;
        return true;
    }
}

template<ui dim> bool md<dim>::rem_typeinteraction(ui type1,ui type2)
{
    pair<ui,ui> id=network.hash(type1,type2);
    if(network.lookup.find(id)!=network.lookup.end())
    {
        ui pos=network.lookup[id];
        network.library[pos]=network.library.back();
        network.backdoor[pos]=network.backdoor.back();
        network.lookup[network.backdoor[pos]]=pos;
        network.library.pop_back();
        network.backdoor.pop_back();
        network.lookup.erase(id);
        return true;
    }
    else return false;
}

template<ui dim> void md<dim>::set_rco(ldf rco)
{
    network.rco=rco;
    network.rcosq=pow(rco,2);
}

template<ui dim> void md<dim>::set_ssz(ldf ssz)
{
    network.ssz=ssz;
    network.sszsq=pow(ssz,2);
}

template<ui dim> ldf md<dim>::dap(ui i,ldf ad)
{
    ldf d;
    switch(simbox.bcond[i])
    {
        case BCOND::PERIODIC: d=fabs(ad)<0.5*simbox.L[i]?ad:ad-fabs(ad+0.5*simbox.L[i])+fabs(ad-0.5*simbox.L[i]); break;
        default: d=ad; break;
    }
    return d;
}

template<ui dim> ldf md<dim>::distsq(ui p1,ui p2)
{
    ldf retval=0.0;
    for(ui i=0;i<dim;i++) retval+=pow(dd(i,p1,p2),2);
    return retval;
}

template<ui dim> ldf md<dim>::dd(ui i,ui p2,ui p1)
{
    ldf ad=particles[p2].x[i]-particles[p1].x[i],d=dap(i,ad);
    return d;
}

template<ui dim> void md<dim>::thread_clear_forces(ui i)
{
    for(ui d=0;d<dim;d++) particles[i].F[d]=0.0;
}

template<ui dim> void md<dim>::thread_calc_forces(ui i)
{
    for(ui j=network.skins[i].size()-1;j<numeric_limits<ui>::max();j--) if(i>network.skins[i][j].neighbor)
    {
        const ldf rsq=distsq(i,network.skins[i][j].neighbor);
        if(!network.update or (network.update and rsq<network.rcosq))
        {
            const ldf r=sqrt(rsq);
            const ldf dVdr=v.dr(network.library[network.skins[i][j].interaction].potential,r,&network.library[network.skins[i][j].interaction].parameters);
            for(ui d=0;d<dim;d++)
            {
                const ldf delta=dd(d,i,network.skins[i][j].neighbor);
                const ldf F=-delta*dVdr/r;
                #ifdef THREADS
                lock_guard<mutex> freeze(parallel.lock);
                particles[i].F[d]-=F;
                particles[network.skins[i][j].neighbor].F[d]+=F;
                #elif OPENMP
                #pragma omp atomic
                particles[i].F[d]-=F;
                #pragma omp atomic
                particles[network.skins[i][j].neighbor].F[d]+=F;
                #else
                particles[i].F[d]-=F;
                particles[network.skins[i][j].neighbor].F[d]+=F;
                #endif
            }
        }
    }
}

template<ui dim> void md<dim>::thread_index_stick(ui i)
{
    for(ui d=0;d<dim;d++) particles[i].xsk[d]=particles[i].x[d];
}

template<ui dim> void md<dim>::index()
{
    #ifdef THREADS
    for(ui t=0;t<parallel.nothreads;t++) parallel.block[t]=thread([=](ui t){for(ui i=t;i<N;i+=parallel.nothreads) thread_index_stick(i);},t);
    for(ui t=0;t<parallel.nothreads;t++) parallel.block[t].join();
    #elif OPENMP
    #pragma omp parallel for
    for(ui i=0;i<N;i++) thread_index_stick(i);;
    #else
    for(ui i=0;i<N;i++) thread_index_stick(i);
    #endif
    switch (indexdata.method)
    {
        case INDEX::BRUTE_FORCE:
			bruteforce();
        break;
		default:
			cell();
        break;
	}
}

template<ui dim> bool md<dim>::test_index()
{
    for(ui i=0;i<N;i++)
    {
        ldf test=0.0;
        for(ui d=0;d<dim;d++) test+=pow(dap(d,particles[i].xsk[d]-particles[i].x[d]),2);
        if(test<pow(network.ssz-network.rco,2)) return true;
    }
    return true;
}

template<ui dim> void md<dim>::calc_forces()
{
    #ifdef THREADS
    for(ui t=0;t<parallel.nothreads;t++) parallel.block[t]=thread([=](ui t){for(ui i=t;i<N;i+=parallel.nothreads) thread_clear_forces(i);},t);
    for(ui t=0;t<parallel.nothreads;t++) parallel.block[t].join();
    for(ui t=0;t<parallel.nothreads;t++) parallel.block[t]=thread([=](ui t){for(ui i=t;i<N;i+=parallel.nothreads) thread_calc_forces(i);},t);
    for(ui t=0;t<parallel.nothreads;t++) parallel.block[t].join();
    #elif OPENMP
    #pragma omp parallel for
    for(ui i=0;i<N;i++) thread_clear_forces(i);
    #pragma omp parallel for
    for(ui i=0;i<N;i++) thread_calc_forces(i);
    #else
    for(ui i=0;i<N;i++) thread_clear_forces(i);
    for(ui i=0;i<N;i++) thread_calc_forces(i);
    #endif
}

template<ui dim> void md<dim>::recalc_forces()
{
    #ifdef THREADS
    for(ui t=0;t<parallel.nothreads;t++) parallel.block[t]=thread([=](ui t){for(ui i=t;i<N;i+=parallel.nothreads) thread_calc_forces(i);},t);
    for(ui t=0;t<parallel.nothreads;t++) parallel.block[t].join();
    #elif OPENMP
    #pragma omp parallel for
    for(ui i=0;i<N;i++) thread_calc_forces(i);
    #else
    for(ui i=0;i<N;i++) thread_calc_forces(i);
    #endif
}

template<ui dim> void md<dim>::thread_periodicity(ui i)
{
    for(ui d=0;d<dim;d++) switch(simbox.bcond[d])
    {
        case BCOND::PERIODIC:
            particles[i].x[d]-=simbox.L[d]*round(particles[i].x[d]/simbox.L[d]);
        break;
        case BCOND::HARD:
            particles[i].x[d]=simbox.L[d]*(fabs(particles[i].x[d]/simbox.L[d]+0.5-2.0*floor(particles[i].x[d]/(2.0*simbox.L[d])+0.75))-0.5);
            particles[i].dx[d]*=-1.0;
        break;
    }
}

template<ui dim> void md<dim>::thread_seuler(ui i)
{
    const ldf o=integrator.h/particles[i].m;
    for(ui d=0;d<dim;d++)
    {
        particles[i].dx[d]-=o*particles[i].F[d];
        particles[i].xp[d]=particles[i].x[d];
        particles[i].x[d]+=integrator.h*particles[i].dx[d];
    }
}

template<ui dim> void md<dim>::thread_vverlet_x(ui i)
{
    for(ui d=0;d<dim;d++)
    {
        particles[i].xp[d]=particles[i].x[d];
        particles[i].x[d]+=integrator.h*particles[i].dx[d]-0.5*pow(integrator.h,2)*particles[i].F[d]/particles[i].m;
    }
}

template<ui dim> void md<dim>::thread_vverlet_dx(ui i)
{
    for(ui d=0;d<dim;d++) particles[i].dx[d]-=0.5*integrator.h*particles[i].F[d]/particles[i].m;
}

template<ui dim> void md<dim>::integrate()
{
    switch(integrator.method)
    {
        case INTEGRATOR::VVERLET:
            #ifdef THREADS
            for(ui t=0;t<parallel.nothreads;t++) parallel.block[t]=thread([=](ui t){for(ui i=t;i<N;i+=parallel.nothreads) if(!particles[i].fix) thread_vverlet_x(i);},t);
            for(ui t=0;t<parallel.nothreads;t++) parallel.block[t].join();
            for(ui t=0;t<parallel.nothreads;t++) parallel.block[t]=thread([=](ui t){for(ui i=t;i<N;i+=parallel.nothreads) thread_calc_forces(i);},t);
            for(ui t=0;t<parallel.nothreads;t++) parallel.block[t].join();
            for(ui t=0;t<parallel.nothreads;t++) parallel.block[t]=thread([=](ui t){for(ui i=t;i<N;i+=parallel.nothreads) if(!particles[i].fix) thread_vverlet_dx(i);},t);
            for(ui t=0;t<parallel.nothreads;t++) parallel.block[t].join();
            #elif OPENMP
            #pragma omp parallel for
            for(ui i=0;i<N;i++) if(!particles[i].fix) thread_vverlet_x(i);
            #pragma omp parallel for
            for(ui i=0;i<N;i++) if(!particles[i].fix) thread_calc_forces(i);
            #pragma omp parallel for
            for(ui i=0;i<N;i++) if(!particles[i].fix) thread_vverlet_dx(i);
            #else
            for(ui i=0;i<N;i++) if(!particles[i].fix) thread_vverlet_x(i);
            recalc_forces();
            for(ui i=0;i<N;i++) if(!particles[i].fix) thread_vverlet_dx(i);
            #endif
        break;
        default:
            #ifdef THREADS
            for(ui t=0;t<parallel.nothreads;t++) parallel.block[t]=thread([=](ui t){for(ui i=t;i<N;i+=parallel.nothreads) if(!particles[i].fix) thread_seuler(i);},t);
            for(ui t=0;t<parallel.nothreads;t++) parallel.block[t].join();
            #elif OPENMP
            #pragma omp parallel for
            for(ui i=0;i<N;i++) if(!particles[i].fix) thread_seuler(i);
            #else
            for(ui i=0;i<N;i++) if(!particles[i].fix) thread_seuler(i);
            #endif
        break;
    }
    #ifdef THREADS
    for(ui t=0;t<parallel.nothreads;t++) parallel.block[t]=thread([=](ui t){for(ui i=t;i<N;i+=parallel.nothreads) if(!particles[i].fix) thread_periodicity(i);},t);
    for(ui t=0;t<parallel.nothreads;t++) parallel.block[t].join();
    #elif OPENMP
    #pragma omp parallel for
    for(ui i=0;i<N;i++) if(!particles[i].fix) thread_periodicity(i);
    #else
    for(ui i=0;i<N;i++) if(!particles[i].fix) thread_periodicity(i);
    #endif
}

template<ui dim> void md<dim>::timestep()
{
    if(network.update and test_index()) index();
    calc_forces();
    integrate();
}

template<ui dim> void md<dim>::timesteps(ui k)
{
    for(ui i=0;i<k;i++) timestep();
}

template<ui dim> ldf md<dim>::thread_H(ui i)
{
    return T(i)+V(i);
}

template<ui dim> ldf md<dim>::thread_T(ui i)
{
    ldf retval=0.0;
    for(ui d=0;d<dim;d++) retval+=pow(particles[i].dx[d],2);
    return 0.5*particles[i].m*retval;
}

template<ui dim> ldf md<dim>::thread_V(ui i)
{
    ldf retval=0.0;
    for(ui j=network.skins[i].size()-1;j<numeric_limits<ui>::max();j--)
    {
        const ldf rsq=distsq(i,network.skins[i][j].neighbor);
        if(rsq<network.rcosq)
        {
            const ldf r=sqrt(r);
            for(ui d=0;d<dim;d++) retval+=(v(network.library[network.skins[i][j].interaction].potential,r,rsq,&network.library[network.skins[i][j].interaction].parameters)-network.library[network.skins[i][j].interaction].vco);
        }
    }
    return retval;
}

template<ui dim> ldf md<dim>::H()
{
    return T()+V();
}

template<ui dim> ldf md<dim>::T()
{
    ldf retval=0.0;
    for(ui i=0;i<N;i++) retval+=thread_T(i);
    return retval/N;
}

template<ui dim> ldf md<dim>::V()
{
    ldf retval=0.0;
    for(ui i=0;i<N;i++) retval+=thread_V(i);
    return retval/N;
}

template<ui dim> void md<dim>::add_particle(ldf mass,ui ptype,bool fixed)
{
    N++;
    particles.push_back(particle<dim>(mass,ptype,fixed));
    network.skins.resize(N);
    index();
}

template<ui dim> void md<dim>::rem_particle(ui particlenr)
{
    N--;
    particles.erase(particlenr);
    network.skins.erase(network.skins.begin()+particlenr);
    index();
}

template<ui dim> void md<dim>::clear()
{
    N=0;
    particles.clear();
    network.skins.clear();
    network.library.clear();
    network.backdoor.clear();
    network.lookup.clear();
}

template<ui dim> void md<dim>::import_pos(...)
{
    ldf *x;
    va_list argv;
    va_start(argv,dim);
    for(ui d=0;d<dim;d++)
    {
        x=va_arg(argv,ldf *);
        #ifdef THREADS
        for(ui t=0;t<parallel.nothreads;t++) parallel.block[t]=thread([=](ui t){for(ui i=t;i<N;i+=parallel.nothreads) particles[i].x[d]=x[i];},t);
        for(ui t=0;t<parallel.nothreads;t++) parallel.block[t].join();
        #elif OPENMP
        #pragma omp parallel for
        for(ui i=0;i<N;i++) particles[i].x[d]=x[i];
        #else
        for(ui i=0;i<N;i++) particles[i].x[d]=x[i];
        #endif
    }
    va_end(argv);
}

template<ui dim> void md<dim>::import_vel(...)
{
    ldf *dx;
    va_list argv;
    va_start(argv,dim);
    for(ui d=0;d<dim;d++)
    {
        dx=va_arg(argv,ldf *);
        #ifdef THREADS
        for(ui t=0;t<parallel.nothreads;t++) parallel.block[t]=thread([=](ui t){for(ui i=t;i<N;i+=parallel.nothreads) particles[i].dx[d]=dx[i];},t);
        for(ui t=0;t<parallel.nothreads;t++) parallel.block[t].join();
        #elif OPENMP
        #pragma omp parallel for
        for(ui i=0;i<N;i++) particles[i].dx[d]=dx[i];
        #else
        for(ui i=0;i<N;i++) particles[i].dx[d]=dx[i];
        #endif
    }
    va_end(argv);
}

template<ui dim> void md<dim>::import_force(...)
{
    ldf *F;
    va_list argv;
    va_start(argv,dim);
    for(ui d=0;d<dim;d++)
    {
        #ifdef THREADS
        for(ui t=0;t<parallel.nothreads;t++) parallel.block[t]=thread([=](ui t){for(ui i=t;i<N;i+=parallel.nothreads) particles[i].F[d]=F[i];},t);
        for(ui t=0;t<parallel.nothreads;t++) parallel.block[t].join();
        #elif OPENMP
        #pragma omp parallel for
        for(ui i=0;i<N;i++) particles[i].F[d]=F[i];
        #else
        for(ui i=0;i<N;i++) particles[i].F[d]=F[i];
        #endif
    }
    va_end(argv);
}

template<ui dim> void md<dim>::export_pos(...)
{
    ldf *x;
    va_list argv;
    va_start(argv,dim);
    for(ui d=0;d<dim;d++)
    {
        x=va_arg(argv,ldf *);
        #ifdef THREADS
        for(ui t=0;t<parallel.nothreads;t++) parallel.block[t]=thread([=](ui t){for(ui i=t;i<N;i+=parallel.nothreads) x[i]=particles[i].x[d];},t);
        for(ui t=0;t<parallel.nothreads;t++) parallel.block[t].join();
        #elif OPENMP
        #pragma omp parallel for
        for(ui i=0;i<N;i++) x[i]=particles[i].x[d];
        #else
        for(ui i=0;i<N;i++) x[i]=particles[i].x[d];
        #endif
    }
    va_end(argv);
}

template<ui dim> void md<dim>::export_vel(...)
{
    ldf *dx;
    va_list argv;
    va_start(argv,dim);
    for(ui d=0;d<dim;d++)
    {
        dx=va_arg(argv,ldf *);
        #ifdef THREADS
        for(ui t=0;t<parallel.nothreads;t++) parallel.block[t]=thread([=](ui t){for(ui i=t;i<N;i+=parallel.nothreads) dx[i]=particles[i].dx[d];},t);
        for(ui t=0;t<parallel.nothreads;t++) parallel.block[t].join();
        #elif OPENMP
        #pragma omp parallel for
        for(ui i=0;i<N;i++) dx[i]=particles[i].dx[d];
        #else
        for(ui i=0;i<N;i++) dx[i]=particles[i].dx[d];
        #endif
    }
    va_end(argv);
}

template<ui dim> void md<dim>::export_force(...)
{
    ldf *F;
    va_list argv;
    va_start(argv,dim);
    for(ui d=0;d<dim;d++)
    {
        F=va_arg(argv,ldf *);
        #ifdef THREADS
        for(ui t=0;t<parallel.nothreads;t++) parallel.block[t]=thread([=](ui t){for(ui i=t;i<N;i+=parallel.nothreads) F[i]=particles[i].F[d];},t);
        for(ui t=0;t<parallel.nothreads;t++) parallel.block[t].join();
        #elif OPENMP
        #pragma omp parallel for
        for(ui i=0;i<N;i++) F[i]=particles[i].F[d];
        #else
        for(ui i=0;i<N;i++) F[i]=particles[i].F[d];
        #endif
    }
    va_end(argv);
}
