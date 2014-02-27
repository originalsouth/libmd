#ifndef libmd_h
#include "../libmd.h"
#endif

template<ui dim> md<dim>::md(ui particlenr)
{
    init(particlenr);
}

template<ui dim> void md<dim>::init(ui particlenr)
{
    N=particlenr;
    DEBUG_2("creating md<%u> with %u particles",dim,N);
    if(N)
    {
        particles.resize(N);
        network.skins.resize(N);
        network.forces.resize(N);
        network.spid.resize(N);
        for(ui i=0;i<N;i++) network.usedtypes[0].insert(i);
        for(ui i=0;i<N;i++) network.spid[i]=std::numeric_limits<ui>::max();
    }
    avars.export_force_calc=true;
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

template<ui dim> ui md<dim>::add_sp_interaction(ui spt,ui p1,ui p2,ui interaction)
{
    if(spt<network.sptypes.size())
    {
        pair<ui,ui> id=network.hash(p1,p2);
        network.sptypes[spt].splookup[id]=interaction;
        return spt;
    }
    else
    {
        spt=network.sptypes.size();
        superparticletype sptype;
        pair<ui,ui> id=network.hash(p1,p2);
        sptype.splookup[id]=interaction;
        network.sptypes.push_back(sptype);
        return spt;
    }
}

template<ui dim> bool md<dim>::mod_sp_interaction(ui spt,ui p1,ui p2,ui interaction)
{
    if(spt<network.sptypes.size())
    {
        add_sp_interaction(spt,p1,p2,interaction);
        return true;
    }
    else return false;
}

template<ui dim> bool md<dim>::rem_sp_interaction(ui spt,ui p1,ui p2)
{
    if(spt<network.sptypes.size())
    {
        pair<ui,ui> id=network.hash(p1,p2);
        if(network.sptypes[spt].splookup.find(id)!=network.sptypes[spt].splookup.end())
        {
            network.sptypes[spt].splookup.erase(id);
            return true;
        }
        else return false;
    }
    else return false;
}

template<ui dim> bool md<dim>::rem_sp_interaction(ui spt)
{
    if(spt<network.sptypes.size())
    {
        ui spn=network.sptypes.size()-1;
        for(ui i=network.superparticles.size();i<numeric_limits<ui>::max();i--)
        {
            if(network.superparticles[i].sptype==spt) network.superparticles[i].sptype=numeric_limits<ui>::max();
            if(network.superparticles[i].sptype==spn) network.superparticles[i].sptype=spt;
        }
        iter_swap(network.superparticles.begin()+spt,network.superparticles.end());
        network.sptypes.pop_back();
        return true;
    }
    else return false;
}

template<ui dim> ui md<dim>::add_forcetype(ui force,vector<ui> *noparticles,vector<ldf> *parameters)
{
    forcetype temp(force,noparticles,parameters);
    network.forcelibrary.push_back(temp);
    return network.forcelibrary.size()-1;
}

template<ui dim> bool md<dim>::mod_forcetype(ui notype,ui force,vector<ui> *noparticles,vector<ldf> *parameters)
{
    if(notype<network.forcelibrary.size())
    {
        forcetype temp(force,noparticles,parameters);
        network.forcelibrary[notype]=temp;
        return true;
    }
    else return false;
}

template<ui dim> bool md<dim>::rem_forcetype(ui notype)
{
    ui pos=network.forcelibrary.size();
    if(notype<pos)
    {
        if(notype==pos-1)
        {
             network.forcelibrary.erase(network.forcelibrary.begin()+notype);
             for(ui i=0;i<N;i++) for(ui j=network.forces[i].size()-1;j<numeric_limits<ui>::max();j--) if(network.forces[i][j]==notype) network.forces[i].erase(network.forces[i].begin()+j);
        }
        else
        {
            network.forcelibrary[notype]=network.forcelibrary[pos-1];
            network.forcelibrary.erase(network.forcelibrary.begin()+pos-1);
            for(ui i=0;i<N;i++) for(ui j=network.forces[i].size()-1;j<numeric_limits<ui>::max();j--)
            {
                if(network.forces[i][j]==notype) network.forces[i].erase(network.forces[i].begin()+j);
                if(network.forces[i][j]==pos-1) network.forces[i][j]=notype;
            }
        }
    }
    else return false;
}

template<ui dim> void md<dim>::assign_forcetype(ui particlenr,ui ftype)
{
    network.forces[particlenr].push_back(ftype);
}

template<ui dim> void md<dim>::assign_all_forcetype(ui ftype)
{
    for(ui i=0;i<N;i++) network.forces[i].push_back(ftype);
}

template<ui dim> void md<dim>::unassign_forcetype(ui particlenr,ui ftype)
{
    for(ui i=network.forces[particlenr].size()-1;i<numeric_limits<ui>::max();i--) if(network.forces[particlenr][i]==ftype)
    {
        network.forces[particlenr].erase(network.forces[particlenr].begin()+i);
        break;
    }
}

template<ui dim> void md<dim>::unassign_all_forcetype(ui ftype)
{
    for(ui j=0;j<N;j++) for(ui i=network.forces[j].size()-1;i<numeric_limits<ui>::max();i++) if(network.forces[j][i]==ftype)
    {
        network.forces[j].erase(network.forces[j].begin()+i);
        break;
    }
}

template<ui dim> void md<dim>::clear_all_assigned_forcetype()
{
    for(ui i=0;i<N;i++) network.forces[i].clear();
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
    if (simbox.boxShear) {
        // use box matrix to calculate distances
        ldf s;
        for (ui j=0;j<dim;j++) {
            s=0;
            //~ printf("\t %d %1.8Lf\n",j,s);
            for (ui k=0;k<dim;k++) {
                //~ printf("%1.8Lf %1.8Lf %1.8Lf %1.8Lf %1.8Lf\n",s,simbox.LshearInv[j][k],particles[p2].x[k],particles[p1].x[k],simbox.LshearInv[j][k]*(particles[p2].x[k]-particles[p1].x[k]));
                s += simbox.LshearInv[j][k]*(particles[p2].x[k]-particles[p1].x[k]);
            }
            if (simbox.bcond[j] == BCOND::PERIODIC || simbox.bcond[j] == BCOND::BOXSHEAR)
                s=fabs(s)<0.5?s:s-fabs(s+0.5)+fabs(s-0.5);
            d += simbox.Lshear[i][j]*s;
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
    if(network.forcelibrary.size() and network.forces[i].size()) for(ui j=network.forces[i].size()-1;j<numeric_limits<ui>::max();j--)
    {
        ui ftype=network.forces[i][j];
        f(network.forcelibrary[ftype].externalforce,&particles[i],nullptr,&network.forcelibrary[ftype].parameters);
    }
}

template<ui dim> void md<dim>::set_index_method(ui method)
{
    indexdata.method=method;
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
    switch(indexdata.method)
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
    avars.export_force_calc=false;
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

template<ui dim> void md<dim>::thread_periodicity_periodic(ui d,ui i)
{
    ldf dx=simbox.L[d]*round(particles[i].x[d]/simbox.L[d]);
    particles[i].xp[d]-=dx;
    particles[i].x[d]-=dx;
}

template<ui dim> void md<dim>::thread_periodicity_boxshear(ui d,ui i) //FIXME: Jayson
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
    ldf xnew=simbox.L[d]*(fabs(particles[i].x[d]/simbox.L[d]+0.5-2.0*floor(particles[i].x[d]/(2.0*simbox.L[d])+0.75))-0.5);
    ldf sign=(((int)round(particles[i].x[d]/simbox.L[d]))&1?-1.0:1.0);
    particles[i].xp[d]+=sign*(xnew-particles[i].x[d]);
    particles[i].x[d]=xnew;
    particles[i].dx[d]*=sign;
}

template<ui dim> void md<dim>::periodicity()
{
    if(simbox.bcond) for(ui d=0;d<dim;d++) switch(simbox.bcond[d])
    {
        case BCOND::PERIODIC:
        {
            #ifdef THREADS
            for(ui t=0;t<parallel.nothreads;t++) parallel.block[t]=thread([=](ui t){for(ui i=t;i<N;i+=parallel.nothreads) if(!particles[i].fix) thread_periodicity_periodic(d,i);},t);
            for(ui t=0;t<parallel.nothreads;t++) parallel.block[t].join();
            #elif OPENMP
            #pragma omp parallel for
            for(ui i=0;i<N;i++) if(!particles[i].fix) thread_periodicity_periodic(d,i);
            #else
            for(ui i=0;i<N;i++) if(!particles[i].fix) thread_periodicity_periodic(d,i);
            #endif
        }
        break;
        case BCOND::BOXSHEAR:
        {
            #ifdef THREADS
            for(ui t=0;t<parallel.nothreads;t++) parallel.block[t]=thread([=](ui t){for(ui i=t;i<N;i+=parallel.nothreads) if(!particles[i].fix) thread_periodicity_boxshear(d,i);},t);
            for(ui t=0;t<parallel.nothreads;t++) parallel.block[t].join();
            #elif OPENMP
            #pragma omp parallel for
            for(ui i=0;i<N;i++) if(!particles[i].fix) thread_periodicity_boxshear(d,i);
            #else
            for(ui i=0;i<N;i++) if(!particles[i].fix) thread_periodicity_boxshear(d,i);
            #endif
        }
        break;
        case BCOND::HARD:
        {
            #ifdef THREADS
            for(ui t=0;t<parallel.nothreads;t++) parallel.block[t]=thread([=](ui t){for(ui i=t;i<N;i+=parallel.nothreads) if(!particles[i].fix) thread_periodicity_hard(d,i);},t);
            for(ui t=0;t<parallel.nothreads;t++) parallel.block[t].join();
            #elif OPENMP
            #pragma omp parallel for
            for(ui i=0;i<N;i++) if(!particles[i].fix) thread_periodicity_hard(d,i);
            #else
            for(ui i=0;i<N;i++) if(!particles[i].fix) thread_periodicity_hard(d,i);
            #endif
        }
        break;
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
    avars.export_force_calc=true;
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
    periodicity();
}

template<ui dim> void md<dim>::update_boundaries()
{   
    // update box matrix for shear
    for(ui j=0;j<dim;j++) for (ui k=0; k<dim; k++)
    {
        simbox.Lshear[j][k] += simbox.vshear[j][k]*integrator.h;
        // shift by appropriate box lengths so that the off-diagonal entries are in the range -L[j][j]/2 to L[j][j]/2 consistent with the positions
        if (j != k)
        {
            while(simbox.Lshear[j][k]>simbox.Lshear[j][j]/2.) simbox.Lshear[j][k]-=simbox.Lshear[j][j];
            while(simbox.Lshear[j][k]<-simbox.Lshear[j][j]/2.) simbox.Lshear[j][k]+=simbox.Lshear[j][j];
        }
    }
    simbox.invert_box();
}

template<ui dim> void md<dim>::timestep()
{
    if(network.update and test_index())
    {
        DEBUG_2("regenerating skinlist");
        index();
    }
    if(simbox.boxShear) update_boundaries();
    calc_forces();
    integrate();
}

template<ui dim> void md<dim>::timesteps(ui k)
{
    for(ui i=0;i<k;i++) timestep();
}

template<ui dim> void md<dim>::set_damping(ldf coefficient)
{
    vector<ldf> parameters(1,coefficient);
    avars.noftypedamping=add_forcetype(EXTFORCE_DAMPING,nullptr,&parameters);
    assign_all_forcetype(avars.noftypedamping);
}

template<ui dim> void md<dim>::unset_damping()
{
    rem_forcetype(avars.noftypedamping);
}

template<ui dim> ldf md<dim>::thread_H(ui i)
{
    return thread_T(i)+thread_V(i);
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
            const ldf r=sqrt(rsq);
            for(ui d=0;d<dim;d++) retval+=(v(network.library[network.skins[i][j].interaction].potential,r,&network.library[network.skins[i][j].interaction].parameters)-network.library[network.skins[i][j].interaction].vco);
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

template<ui dim> ui md<dim>::add_particle(ldf mass,ui ptype,bool fixed)
{
    N++;
    particles.push_back(particle<dim>(mass,ptype,fixed));
    network.spid.push_back(numeric_limits<ui>::max());
    network.skins.resize(N);
    network.forces.resize(N);
    network.usedtypes[ptype].insert(N-1);
    index();
    return N-1;
}

template<ui dim> ui md<dim>::add_particle(ldf x[dim],ldf mass,ui ptype,bool fixed)
{
    ui i=add_particle(mass,ptype,fixed);
    DEBUG_2("Created particle #%u",i);
    for(ui d=0;d<dim;d++)
    {
        particles[i].x[d]=x[d];
        particles[i].dx[d]=0.0;
    }
    return i;
}

template<ui dim> ui md<dim>::add_particle(ldf x[dim],ldf dx[dim],ldf mass,ui ptype,bool fixed)
{
    ui i=add_particle(mass,ptype,fixed);
    DEBUG_2("Created particle #%u",i);
    for(ui d=0;d<dim;d++)
    {
        particles[i].x[d]=x[d];
        particles[i].dx[d]=dx[d];
    }
    return i;
}

template<ui dim> void md<dim>::rem_particle(ui i)
{
    if(network.spid[i]<N) sp_dispose(network.spid[i],i);
    N--;
    // store particle types of deleted particle and particle that will replace it
    ui deleted_ptype = particles[i].type;
    ui last_ptype = particles.rbegin()->type;
    // swap particle to delete with last particle, to prevent changing index of all particles after particlenr, and then delete it
    std::iter_swap(particles.begin()+i, particles.rbegin());
    particles.pop_back();
    // update usedtypes dictionary
    network.usedtypes[deleted_ptype].erase(i);
    network.usedtypes[last_ptype].erase(N);
    network.usedtypes[last_ptype].insert(i);
    // update the network  TODO: Benny and Thomas -- please check that no other data structures need updating
    std::iter_swap(network.skins.begin()+i, network.skins.rbegin());
    network.skins.pop_back();
    std::iter_swap(network.forces.begin()+i, network.forces.rbegin());
    network.forces.pop_back();
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
    network.spid.clear();
    network.superparticles.clear();
    network.sptypes.clear();
    network.usedtypes.clear();
    network.forcelibrary.clear();
    network.forces.clear();
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

template<ui dim> void md<dim>::import_pos(ui i,ldf x)
{
    ui d=vvars[0];
    particles[i].x[d]=x;
}

template<ui dim> template<typename...arg> void md<dim>::import_pos(ui i,ldf x,arg...argv)
{
    ui d=vvars[0];
    particles[i].x[d]=x;
    import_pos(i,argv...);
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

template<ui dim> void md<dim>::import_vel(ui i,ldf dx)
{
    ui d=vvars[1];
    particles[i].dx[d]=dx;
}

template<ui dim> template<typename...arg> void md<dim>::import_vel(ui i,ldf dx,arg...argv)
{
    ui d=vvars[1];
    particles[i].dx[d]=dx;
    import_vel(i,argv...);
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

template<ui dim> void md<dim>::import_force(ui i,ldf F)
{
    ui d=vvars[2];
    particles[i].F[d]=F;
}

template<ui dim> template<typename...arg> void md<dim>::import_force(ui i,ldf F,arg...argv)
{
    ui d=vvars[2];
    particles[i].F[d]=F;
    import_force(i,argv...);
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
    export_pos(argv...);
}

template<ui dim> void md<dim>::export_pos(ui i,ldf &x)
{
    ui d=vvars[3];
    x=particles[i].x[d];
}

template<ui dim> template<typename...arg> void md<dim>::export_pos(ui i,ldf &x,arg...argv)
{
    ui d=vvars[3];
    x=particles[i].x[d];
    export_pos(i,argv...);
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
    export_vel(argv...);
}

template<ui dim> void md<dim>::export_vel(ui i,ldf &dx)
{
    ui d=vvars[4];
    dx=particles[i].dx[d];
}

template<ui dim> template<typename...arg> void md<dim>::export_vel(ui i,ldf &dx,arg...argv)
{
    ui d=vvars[4];
    dx=particles[i].dx[d];
    export_vel(i,argv...);
}

template<ui dim> void md<dim>::export_force(ldf *F)
{
    ui d=vvars[5];
    if(avars.export_force_calc) calc_forces();
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
    if(avars.export_force_calc) calc_forces();
    #ifdef THREADS
    for(ui t=0;t<parallel.nothreads;t++) parallel.block[t]=thread([=](ui t){for(ui i=t;i<N;i+=parallel.nothreads) F[i]=particles[i].F[d];},t);
    for(ui t=0;t<parallel.nothreads;t++) parallel.block[t].join();
    #elif OPENMP
    #pragma omp parallel for
    for(ui i=0;i<N;i++) F[i]=particles[i].F[d];
    #else
    for(ui i=0;i<N;i++) F[i]=particles[i].F[d];
    #endif
    export_force(argv...);
}

template<ui dim> void md<dim>::export_force(ui i,ldf &F)
{
    ui d=vvars[5];
    if(avars.export_force_calc) calc_forces();
    F=particles[i].F[d];
}

template<ui dim> template<typename...arg> void md<dim>::export_force(ui i,ldf &F,arg...argv)
{
    ui d=vvars[5];
    if(!d) calc_forces();
    F=particles[i].F[d];
    export_force(i,argv...);
}

template<ui dim> ldf md<dim>::direct_readout(ui i,uc type)
{
    ui d=vvars[6];
    switch(type)
    {
        case 'v': return particles[i].dx[d]; break;
        case 'F':
        {
            if(avars.export_force_calc) calc_forces();
            return particles[i].F[d];
        }
        break;
        default: return particles[i].x[d]; break;
    }
}

template<ui dim> ldf md<dim>::direct_readout(ui d,ui i,uc type)
{
    switch(type)
    {
        case 'v': return particles[i].dx[d]; break;
        case 'F':
        {
            if(avars.export_force_calc) calc_forces();
            return particles[i].F[d];
        }
        break;
        default: return particles[i].x[d]; break;
    }
}

template<ui dim> void md<dim>::fix_particle(ui i,bool fix)
{
    DEBUG_2("Fixing(%d) particle %u.",fix,i);
    particles[i].fix=fix;
}

template<ui dim> void md<dim>::fix_particles(ui spi,bool fix)
{
    DEBUG_2("Fixing(%d) super particle particle %u.",fix,spi);
    for(auto it=network.superparticles[spi].particles.begin();it!=network.superparticles[spi].particles.end();it++) particles[it->first].fix=fix;
}

template<ui dim> void md<dim>::translate_particle(ui i,ldf x[dim])
{
    DEBUG_2("Translating particle particle %u.",i);
    for(ui d=0;d<dim;d++) particles[i].x[d]+=x[d];
}

template<ui dim> void md<dim>::translate_particles(ui spi,ldf x[dim])
{
    DEBUG_2("Translating super particle particle %u.",spi);
    for(auto it=network.superparticles[spi].particles.begin();it!=network.superparticles[spi].particles.end();it++) for(ui d=0;d<dim;d++) particles[it->first].x[d]+=x[d];
}

template<ui dim> void md<dim>::drift_particle(ui i,ldf dx[dim])
{
    DEBUG_2("Drifting particle particle %u.",i);
    for(ui d=0;d<dim;d++) particles[i].dx[d]+=dx[d];
}

template<ui dim> void md<dim>::drift_particles(ui spi,ldf dx[dim])
{
    DEBUG_2("Drifting Translating super particle particle %u.",spi);
    for(auto it=network.superparticles[spi].particles.begin();it!=network.superparticles[spi].particles.end();it++) for(ui d=0;d<dim;d++) particles[it->first].dx[d]+=dx[d];
}

template<ui dim> void md<dim>::set_position_particles(ui spi,ldf x[dim])
{
    DEBUG_2("Drifting Translating super particle particle %u.",spi);
    ldf dx[dim];
    get_position_particles(spi,dx);
    for(ui d=0;d<dim;d++) dx[d]=x[d]-dx[d];
    for(auto it=network.superparticles[spi].particles.begin();it!=network.superparticles[spi].particles.end();it++) for(ui d=0;d<dim;d++) particles[it->first].x[d]+=dx[d];
}

template<ui dim> void md<dim>::set_velocity_particles(ui spi,ldf dx[dim])
{
    DEBUG_2("Drifting Translating super particle particle %u.",spi);
    for(auto it=network.superparticles[spi].particles.begin();it!=network.superparticles[spi].particles.end();it++) for(ui d=0;d<dim;d++) particles[it->first].dx[d]=dx[d];
}

template<ui dim> void md<dim>::get_position_particles(ui spi,ldf x[dim])
{
    DEBUG_2("Calculating center of mass super particle particle %u.",spi);
    ldf m=0.0;
    for(ui d=0;d<dim;d++) x[d]=0.0;
    for(auto it=network.superparticles[spi].particles.begin();it!=network.superparticles[spi].particles.end();it++) for(ui d=0;d<dim;d++)
    {
        x[d]+=particles[it->first].m*particles[it->first].x[d];
        m+=particles[it->first].m;
    }
    for(ui d=0;d<dim;d++) x[d]/=m;
}

template<ui dim> void md<dim>::get_velocity_particles(ui spi,ldf dx[dim])
{
    DEBUG_2("Drifting Translating super particle particle %u.",spi);
    for(ui d=0;d<dim;d++) dx[d]=0.0;
    for(auto it=network.superparticles[spi].particles.begin();it!=network.superparticles[spi].particles.end();it++) for(ui d=0;d<dim;d++) dx[d]+=particles[it->first].dx[d];
    for(ui d=0;d<dim;d++) dx[d]/=network.superparticles[spi].particles.size();
}

template<ui dim> ui md<dim>::sp_ingest(ui spi,ui i)
{
    if(spi<network.superparticles.size())
    {
        if(network.spid[i]==spi)
        {
            WARNING("Praticle %u is already in super_particle %u.",i,spi);
            return spi;
        }
        network.spid[i]=spi;
        network.superparticles[spi].particles[i]=network.superparticles[spi].particles.size();
    }
    else
    {
        spi=network.superparticles.size();
        network.spid[i]=spi;
        superparticle sp;
        sp.particles[i]=sp.particles.size();
        sp.sptype=numeric_limits<ui>::max();
        network.superparticles.push_back(sp);
    }
    return spi;
}

template<ui dim> ui md<dim>::sp_ingest(ui spi,ui sptype,ui i)
{
    if(spi<network.superparticles.size())
    {
        network.spid[i]=spi;
        network.superparticles[spi].particles[i]=network.superparticles[spi].particles.size();
        network.superparticles[spi].sptype=sptype;
    }
    else
    {
        spi=network.superparticles.size();
        network.spid[i]=spi;
        superparticle sp;
        sp.particles[i]=sp.particles.size();
        sp.sptype=sptype;
        network.superparticles.push_back(sp);
    }
    return spi;
}

template<ui dim> void md<dim>::sp_dispose(ui spi)
{
    if(spi<network.superparticles.size())
    {
        ui spn=network.superparticles.size()-1;
        if(spi<spn)
        {
            for(auto it=network.superparticles[spn].particles.begin();it!=network.superparticles[spn].particles.end();it++) network.spid[it->first]=spi;
            iter_swap(network.superparticles.begin()+spi,network.superparticles.end());
        }
        for(auto it=network.superparticles[spn].particles.begin();it!=network.superparticles[spn].particles.end();it++) network.spid[it->first]=numeric_limits<ui>::max();
        network.superparticles.pop_back();
    }
}

template<ui dim> void md<dim>::sp_p_dispose(ui i)
{
    if(i<N)
    {
        ui spi=network.spid[i];
        if(spi<network.superparticles.size())
        {
            if(network.superparticles[spi].particles.size()<2) sp_dispose(spi);
            else
            {
                network.spid[i]=numeric_limits<ui>::max();
                network.superparticles[spi].particles.erase(i);
            }
        }
    }
}

template<ui dim> void md<dim>::add_bond(ui p1, ui p2, ui itype, vector<ldf> *params)
{
    /* add a 'bond' i.e. a specific interaction between two particles, of type itype and with parameter params */
    /* NOTE: forces p1 and p2 to have unique particle types. Replicates former interactions experienced between
     * p1 or p2 and other particle types. */
    
    set<ui> partners_of_p1, partners_of_p2;
    
    // assign unique types to points
    ui old_p1type = particles[p1].type;
    ui new_p1type = old_p1type;
    
    if (network.usedtypes[old_p1type].size() > 1) {
        // p1 does not have a unique type; reassign. 
        new_p1type = (network.usedtypes.rbegin()->first + 1); // use one greater than largest used particle type
        set_type(p1, new_p1type);
        
        // keep track of previously defined interactions
        for (map<pair<ui,ui>,ui>::iterator it = network.lookup.begin(); it != network.lookup.end(); it++) {
            pair<ui, ui> ipair = it->first;
            if (ipair.first == old_p1type) partners_of_p1.insert(ipair.second);
            else if (ipair.second == old_p1type) partners_of_p1.insert(ipair.first);
        }
    }
    
    ui old_p2type = particles[p2].type;
    ui new_p2type = old_p2type;
    
    if (network.usedtypes[old_p2type].size() > 1) {
        // p2 does not have a unique type; reassign.
        new_p2type = (network.usedtypes.rbegin()->first + 1); // use one greater than largest used particle type
        set_type(p2, new_p2type);
        
        // keep track of previously defined interactions
        for (map<pair<ui,ui>,ui>::iterator it = network.lookup.begin(); it != network.lookup.end(); it++) {
            pair<ui, ui> ipair = it->first;
            if (ipair.first == old_p2type) partners_of_p2.insert(ipair.second);
            else if (ipair.second == old_p2type) partners_of_p2.insert(ipair.first);
        }
    }
    
    // now add the interaction
    rem_typeinteraction(new_p1type, new_p2type); // removes any old interaction between the unique ids new_p1type and new_p2type
    add_typeinteraction(new_p1type, new_p2type, itype, params);
    
    // loop through previously defined interactions and clone them so that they are preserved
    if (partners_of_p1.size() > 0) {
        for (set<ui>::iterator it = partners_of_p1.begin(); it != partners_of_p1.end(); it++) {
            if (*it != new_p2type) {
                interactiontype old_interaction = network.library[network.lookup[network.hash(old_p1type,*it)]];
                vector<ldf> old_params = old_interaction.parameters;
                add_typeinteraction(new_p1type, *it, old_interaction.potential, &old_params);
            } 
        }
    }
    if (partners_of_p2.size() > 0) {
        for (set<ui>::iterator it = partners_of_p2.begin(); it != partners_of_p2.end(); it++) {
            if (*it != new_p1type) {
                interactiontype old_interaction = network.library[network.lookup[network.hash(old_p2type,*it)]];
                vector<ldf> old_params = network.library[network.lookup[network.hash(old_p2type,*it)]].parameters;
                add_typeinteraction(new_p2type, *it, old_interaction.potential, &old_params);
            } 
        }
    }
    
}

template<ui dim> void md<dim>::add_spring(ui p1, ui p2, ldf springconstant, ldf l0)
{
    /* add a spring between two points with specified springconstant and equilibrium length */
    vector<ldf> params = {springconstant, l0};
    add_bond(p1,p2,POT::POT_HOOKEAN,&params);
}

template<ui dim> bool md<dim>::share_bond(ui p1, ui p2)
{
    /* Check whether particles p1 and p2 share a bond. */
    
    // 1. Do the particles have unique types?
    if (network.usedtypes[particles[p1].type].size() > 1 || network.usedtypes[particles[p2].type].size() > 1) return false;
    
    // 2. Do the unique types have an interaction entry?
    pair<ui,ui> id=network.hash(particles[p1].type,particles[p2].type);
    if(network.lookup.find(id)==network.lookup.end()) return false;
    
    // bond exists. NOTE: Does not take indexing into account. TODO?
    return true;
}

template<ui dim> bool md<dim>::rem_bond(ui p1, ui p2)
{
    /* remove bond-style interaction between particles p1 and p2. does not affect other interactions. */
    if (!share_bond(p1, p2)) return false;
    return rem_typeinteraction(particles[p1].type, particles[p2].type);
}

template<ui dim> bool md<dim>::mod_bond(ui p1, ui p2, ui potential, vector<ldf> *parameters)
{
    /* modify bond-style interaction between particles p1 and p2. does not affect other interactions. */
    if (!share_bond(p1, p2)) return false;
    return mod_typeinteraction(particles[p1].type, particles[p2].type, potential, parameters);
}

template<ui dim> void md<dim>::set_type(ui p, ui newtype)
{
    /* change a particle's type and update  network.usedtypes */
    ui oldtype = particles[p].type;
    if (oldtype != newtype) {
        particles[p].type = newtype;
        network.usedtypes[oldtype].erase(p);
        network.usedtypes[newtype].insert(p);
    }
}


