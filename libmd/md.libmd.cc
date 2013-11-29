template<ui dim> md<dim>::md()
{
    N=0;
    nothreads=1;
}

template<ui dim> md<dim>::md(ui particlenr)
{
    N=particlenr;
    nothreads=1;
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
        case 1: d=fabs(ad)<0.5*simbox.L[i]?ad:ad-fabs(ad+0.5*simbox.L[i])+fabs(ad-0.5*simbox.L[i]); break;
        default: d=ad; break;
    }
    return d;
}

template<ui dim> void md<dim>::thread_clear_forces(ui i)
{
    for(ui d=0;d<dim;d++) particles[i].F[d]=0.0;
}

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
                lock_guard<mutex> freeze(lock);
                particles[i].F[d]+=F;
                particles[network.skins[i][j].neighbor].F[d]-=F;
            }
        }
    }
}

template<ui dim> void md<dim>::index()
{
    for(ui i=0;i<N;i++)
    {
        network.skins[i].clear();
        for(ui j=0;j<N;j++) if(i!=j and distsq(i,j)<network.sszsq)
        {
            const pair<ui,ui> it=network.hash(particles[i].type,particles[j].type);
            if(network.lookup.count(it))
            {
                interactionneighbor in(j,network.lookup[it]);
                network.skins[i].push_back(in);
            }
        }
    }
}

template<ui dim> void md<dim>::calc_forces()
{
    const ui max=nothreads;
    vector<thread> block(max);
    for(ui t=0;t<max;t++) block[t]=thread([=](ui t){for(ui i=t;i<N;i+=max) thread_clear_forces(i);},t);
    for(ui t=0;t<max;t++) block[t].join();
    for(ui t=0;t<max;t++) block[t]=thread([=](ui t){for(ui i=t;i<N;i+=max) thread_calc_forces(i);},t);
    for(ui t=0;t<max;t++) block[t].join();
}

template<ui dim> void md<dim>::recalc_forces()
{
    const ui max=nothreads;
    vector<thread> block(max);
    for(ui t=0;t<max;t++) block[t]=thread([=](ui t){for(ui i=t;i<N;i+=max) thread_calc_forces(i);},t);
    for(ui t=0;t<max;t++) block[t].join();
}

template<ui dim> void md<dim>::thread_periodicity(ui i)
{
    for(ui d=0;d<dim;d++) switch(simbox.bcond[d])
    {
        case 1:
            particles[i].x[d]-=simbox.L[d]*floor((particles[i].x[d]+simbox.L[d]/2.0)/simbox.L[d]);
        break;
        case 2:
            particles[i].x[d]-=simbox.L[d]*round(particles[i].x[d]/simbox.L[d]);
            particles[i].dx[d]*=-1.0;
        break;
    }
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

template<ui dim> void md<dim>::integrate()
{
    const ui max=nothreads;
    vector<thread> block(max);
    switch(integrator.method)
    {
        case 1:
            for(ui t=0;t<max;t++) block[t]=thread([=](ui t){for(ui i=t;i<N;i+=max) if(!particles[i].fix) thread_integrate(i,0);},t);
            for(ui t=0;t<max;t++) block[t].join();
            for(ui t=0;t<max;t++) block[t]=thread([=](ui t){for(ui i=t;i<N;i+=max) thread_calc_forces(i);},t);
            for(ui t=0;t<max;t++) block[t].join();
            for(ui t=0;t<max;t++) block[t]=thread([=](ui t){for(ui i=t;i<N;i+=max) if(!particles[i].fix) thread_integrate(i,1);},t);
            for(ui t=0;t<max;t++) block[t].join();
        break;
        default:
            for(ui t=0;t<max;t++) block[t]=thread([=](ui t){for(ui i=t;i<N;i+=max) if(!particles[i].fix) thread_integrate(i,0);},t);
            for(ui t=0;t<max;t++) block[t].join();
        break;
    }
    for(ui t=0;t<max;t++) block[t]=thread([=](ui t){for(ui i=t;i<N;i+=max) thread_periodicity(i);},t);
    for(ui t=0;t<max;t++) block[t].join();
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

//TODO: Test if this is faster than summing over H(i)
template<ui dim> ldf md<dim>::H()
{
    return T()+V();
}

//TODO: Make parallel launcher
template<ui dim> ldf md<dim>::T()
{
    ldf retval=0.0;
    for(ui i=0;i<N;i++) retval+=T(i);
    return retval/N;
}

//TODO: Make parallel launcher
template<ui dim> ldf md<dim>::V()
{
    ldf retval=0.0;
    for(ui i=0;i<N;i++) retval+=V(i);
    return retval/N;
}
