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
    for (ui i = 0; i < N; i++) network.usedtypes[0].insert(i); // assumes default particle type is 0. fix later?
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

template<ui dim> ldf md<dim>::dd(ui i,ui p1,ui p2) //TODO: fix non-periodic boundary conditions plus shear; fix names
{   
    ldf d=0;
    if (simbox.LeesEdwards) {
        // use box matrix to calculate distances
        ldf sab[dim] = {};
        for (ui j=0;j<dim;j++) {
            sab[j]=0;
            //~ printf("\t %d %1.8Lf\n",j,sab[j]);
            for (ui k=0;k<dim;k++) {
                //~ printf("%1.8Lf %1.8Lf %1.8Lf %1.8Lf %1.8Lf\n",sab[j],simbox.LshearInv[j][k],particles[p2].x[k],particles[p1].x[k],simbox.LshearInv[j][k]*(particles[p2].x[k]-particles[p1].x[k]));
                sab[j] += simbox.LshearInv[j][k]*(particles[p2].x[k]-particles[p1].x[k]);
            }
            sab[j]=fabs(sab[j])<0.5?sab[j]:sab[j]-fabs(sab[j]+0.5)+fabs(sab[j]-0.5);
            d += simbox.Lshear[i][j]*sab[j];
        }
    }
    else {
        ldf ad=particles[p2].x[i]-particles[p1].x[i];
        d=dap(i,ad);
    }
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
                const ldf F_i=dd(d,i,network.skins[i][j].neighbor)*dVdr/r;
                #ifdef THREADS
                lock_guard<mutex> freeze(parallel.lock);
                particles[i].F[d]+=F_i;
                particles[network.skins[i][j].neighbor].F[d]-=F_i;
                #elif OPENMP
                #pragma omp atomic
                particles[i].F[d]+=F_i;
                #pragma omp atomic
                particles[network.skins[i][j].neighbor].F[d]-=F_i;
                #else
                particles[i].F[d]+=F_i;
                particles[network.skins[i][j].neighbor].F[d]-=F_i;
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
    if (simbox.LeesEdwards) {
        // ignore bcond[d] for now; assume all periodic and have entries in Lshear
        // TODO: allow mixed boundary conditions: periodic along directions required by lees-edwards; aperiodic otherwise
        for(ui j=0;j<dim;j++) {
            ldf boundaryCrossing = round(particles[i].x[j]/simbox.L[j]);
            if (fabs(boundaryCrossing) > .1){
                for (ui k=0; k<dim; k++) {
                    particles[i].x[k] -= simbox.Lshear[k][j]*boundaryCrossing;
                    particles[i].dx[k] -= simbox.vshear[k][j]*boundaryCrossing;
                }
            }
        }
    }
    else {
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
}

template<ui dim> void md<dim>::thread_seuler(ui i)
{
    const ldf o=integrator.h/particles[i].m;
    for(ui d=0;d<dim;d++)
    {
        particles[i].dx[d]+=o*particles[i].F[d];
        particles[i].xp[d]=particles[i].x[d];
        particles[i].x[d]+=integrator.h*particles[i].dx[d];
    }
}

template<ui dim> void md<dim>::thread_vverlet_x(ui i)
{
    for(ui d=0;d<dim;d++)
    {
        particles[i].xp[d]=particles[i].x[d];
        particles[i].x[d]+=integrator.h*particles[i].dx[d]+0.5*pow(integrator.h,2)*particles[i].F[d]/particles[i].m;
    }
}

template<ui dim> void md<dim>::thread_vverlet_dx(ui i)
{
    for(ui d=0;d<dim;d++) particles[i].dx[d]+=0.5*integrator.h*particles[i].F[d]/particles[i].m;
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

template<ui dim> void md<dim>::update_boundaries()
{   
    // update box matrix for Lees-Edwards shear
    // TODO: test robustness for shear in multiple directions/boundaries
    for(ui j=0;j<dim;j++) {
        for (ui k=0; k<dim; k++) {
            simbox.Lshear[j][k] += simbox.vshear[j][k]*integrator.h;
            // shift by appropriate box lengths so that the off-diagonal entries are in the range -L[j][j]/2 to L[j][j]/2 consistent with the positions
            if (j != k) {
                while (simbox.Lshear[j][k] > simbox.Lshear[j][j]/2.) simbox.Lshear[j][k] -= simbox.Lshear[j][j];
                while (simbox.Lshear[j][k] <- simbox.Lshear[j][j]/2.) simbox.Lshear[j][k] += simbox.Lshear[j][j];
            }
        }
    }
    simbox.invert_box();
}

template<ui dim> void md<dim>::timestep()
{
    if(network.update and test_index()) index();
    if (simbox.LeesEdwards) update_boundaries();
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
    network.usedtypes[ptype].insert(N-1);
    index();
}

template<ui dim> void md<dim>::rem_particle(ui particlenr)
{
    N--;
    
    // store particle types of deleted particle and particle that will replace it
    ui deleted_ptype = particles[particlenr].type;
    ui last_ptype = particles.rbegin()->type;
    
    // swap particle to delete with last particle, to prevent changing index of all particles after particlenr, and then delete it
    std::iter_swap(particles.begin()+particlenr, particles.rbegin()); 
    particles.pop_back();
    
    // update usedtypes dictionary
    network.usedtypes[deleted_ptype].erase(particlenr);
    network.usedtypes[last_ptype].erase(N);
    network.usedtypes[last_ptype].insert(particlenr);

    // update the network  TODO: benny and thomas -- please check that no other data structures need updating
    std::iter_swap(network.skins.begin()+particlenr, network.skins.rbegin());
    network.skins.pop_back();
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

template<ui dim> void md<dim>::import_pos(ldf *x)
{
    ui d=vvars[0];
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

template<ui dim> template<typename...arg> void md<dim>::import_pos(ldf *x,arg...argv)
{
    ui d=vvars[0];
    #ifdef THREADS
    for(ui t=0;t<parallel.nothreads;t++) parallel.block[t]=thread([=](ui t){for(ui i=t;i<N;i+=parallel.nothreads) particles[i].x[d]=x[i];},t);
    for(ui t=0;t<parallel.nothreads;t++) parallel.block[t].join();
    #elif OPENMP
    #pragma omp parallel for
    for(ui i=0;i<N;i++) particles[i].x[d]=x[i];
    #else
    for(ui i=0;i<N;i++) particles[i].x[d]=x[i];
    #endif
    import_pos(argv...);
}

template<ui dim> void md<dim>::import_vel(ldf *dx)
{
    ui d=vvars[1];
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

template<ui dim> template<typename...arg> void md<dim>::import_vel(ldf *dx,arg...argv)
{
    ui d=vvars[1];
    #ifdef THREADS
    for(ui t=0;t<parallel.nothreads;t++) parallel.block[t]=thread([=](ui t){for(ui i=t;i<N;i+=parallel.nothreads) particles[i].dx[d]=dx[i];},t);
    for(ui t=0;t<parallel.nothreads;t++) parallel.block[t].join();
    #elif OPENMP
    #pragma omp parallel for
    for(ui i=0;i<N;i++) particles[i].dx[d]=dx[i];
    #else
    for(ui i=0;i<N;i++) particles[i].dx[d]=dx[i];
    #endif
    import_vel(argv...);
}

template<ui dim> void md<dim>::import_force(ldf *F)
{
    ui d=vvars[2];
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

template<ui dim> template<typename...arg> void md<dim>::import_force(ldf *F,arg...argv)
{
    ui d=vvars[2];
    #ifdef THREADS
    for(ui t=0;t<parallel.nothreads;t++) parallel.block[t]=thread([=](ui t){for(ui i=t;i<N;i+=parallel.nothreads) particles[i].F[d]=F[i];},t);
    for(ui t=0;t<parallel.nothreads;t++) parallel.block[t].join();
    #elif OPENMP
    #pragma omp parallel for
    for(ui i=0;i<N;i++) particles[i].F[d]=F[i];
    #else
    for(ui i=0;i<N;i++) particles[i].F[d]=F[i];
    #endif
    import_force(argv...);
}

template<ui dim> void md<dim>::export_pos(ldf *x)
{
    ui d=vvars[3];
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

template<ui dim> template<typename...arg> void md<dim>::export_pos(ldf *x,arg...argv)
{
    ui d=vvars[3];
    #ifdef THREADS
    for(ui t=0;t<parallel.nothreads;t++) parallel.block[t]=thread([=](ui t){for(ui i=t;i<N;i+=parallel.nothreads) x[i]=particles[i].x[d];},t);
    for(ui t=0;t<parallel.nothreads;t++) parallel.block[t].join();
    #elif OPENMP
    #pragma omp parallel for
    for(ui i=0;i<N;i++) x[i]=particles[i].x[d];
    #else
    for(ui i=0;i<N;i++) x[i]=particles[i].x[d];
    #endif
    import_pos(argv...);
}

template<ui dim> void md<dim>::export_vel(ldf *dx)
{
    ui d=vvars[4];
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

template<ui dim> template<typename...arg> void md<dim>::export_vel(ldf *dx,arg...argv)
{
    ui d=vvars[4];
    #ifdef THREADS
    for(ui t=0;t<parallel.nothreads;t++) parallel.block[t]=thread([=](ui t){for(ui i=t;i<N;i+=parallel.nothreads) dx[i]=particles[i].dx[d];},t);
    for(ui t=0;t<parallel.nothreads;t++) parallel.block[t].join();
    #elif OPENMP
    #pragma omp parallel for
    for(ui i=0;i<N;i++) dx[i]=particles[i].dx[d];
    #else
    for(ui i=0;i<N;i++) dx[i]=particles[i].dx[d];
    #endif
    import_pos(argv...);
}

template<ui dim> void md<dim>::export_force(ldf *F)
{
    ui d=vvars[5];
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

template<ui dim> template<typename...arg> void md<dim>::export_force(ldf *F,arg...argv)
{
    ui d=vvars[5];
    #ifdef THREADS
    for(ui t=0;t<parallel.nothreads;t++) parallel.block[t]=thread([=](ui t){for(ui i=t;i<N;i+=parallel.nothreads) F[i]=particles[i].F[d];},t);
    for(ui t=0;t<parallel.nothreads;t++) parallel.block[t].join();
    #elif OPENMP
    #pragma omp parallel for
    for(ui i=0;i<N;i++) F[i]=particles[i].F[d];
    #else
    for(ui i=0;i<N;i++) F[i]=particles[i].F[d];
    #endif
    import_pos(argv...);
}

template<ui dim> void md<dim>::add_bond(ui p1, ui p2, ui itype, vector<ldf> *params) {
    //TODO: Reassign original interaction to the new particle types

    /* add a 'bond' i.e. a specific interaction between two particles, of type itype and with parameter params */
    /* NOTE: forces p1 and p2 to have unique particle types. Consequently, any interactions experienced by p1 and p2 that involved shared types with other particles.
     * however, bond-like interactions between p1/p2 and other particles are preserved, because they already provided a unique particle type
     * To maintain bonds on top of previous interactions, consider multiple types for each particle */
    
    // assign unique types to points
    ui old_p1type = particles[p1].type;
    ui new_p1type = old_p1type;
    
    if (network.usedtypes[old_p1type].size() > 1) {
        // p1 does not have a unique type; reassign. this undoes whatever interaction p1 used to have.
        new_p1type = (network.usedtypes.rbegin()->first + 1); // use one greater than largest used particle type
        set_type(p1, new_p1type);
    }
    
    ui old_p2type = particles[p2].type;
    ui new_p2type = old_p2type;
    
    if (network.usedtypes[old_p2type].size() > 1) {
        // p2 does not have a unique type; reassign. this undoes whatever interaction p2 used to have.
        new_p2type = (network.usedtypes.rbegin()->first + 1); // use one greater than largest used particle type
        set_type(p2, new_p2type);
    }
    
    // now add the interaction
    add_typeinteraction(new_p1type, new_p2type, itype, params);
    
    // loop through previous interactions and clone them so that they are preserved
    
}

template<ui dim> void md<dim>::add_spring(ui p1, ui p2, ldf springconstant, ldf l0) {
    /* add a spring between two points with specified springconstant and equilibrium length */
    vector<ldf> params = {springconstant, l0};
    add_bond(p1,p2,POT::POT_HOOKIAN,&params);
}


template<ui dim> void md<dim>::set_type(ui p, ui newtype) {
    /* change a particle's type and update  network.usedtypes */
    ui oldtype = particles[p].type;
    if (oldtype != newtype) {
        particles[p].type = newtype;
        network.usedtypes[oldtype].erase(p);
        network.usedtypes[newtype].insert(p);
    }
}


