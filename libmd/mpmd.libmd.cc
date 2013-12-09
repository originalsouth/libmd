template<ui dim> mpmd<dim>::mpmd():md<dim>()
{

}

template<ui dim> mpmd<dim>::mpmd(ui particlenr):md<dim>(particlenr)
{

}

template<ui dim> ldf mpmd<dim>::embedded_distsq(ui p1,ui p2)
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
    ldf ad=pow(patch.f(particles[p2].x)-patch.f(particles[p1].x),2);
    return retval+ad;
}

template<ui dim> ldf mpmd<dim>::embedded_dd(ui i,ui p1,ui p2)
{
    ldf ad=particles[p2].x[i]-particles[p1].x[i],d;
    switch(simbox.bcond[i])
    {
        case 1: d=fabs(ad)<0.5*simbox.L[i]?ad:ad-fabs(ad+0.5*simbox.L[i])+fabs(ad-0.5*simbox.L[i]); break;
        default: d=ad; break;
    }
    ad=patch.f(particles[p2].x)-patch.f(particles[p1].x);
    return d+ad;
}

template<ui dim> ldf mpmd<dim>::gimmel(ui i,ui s)
{
    ldf retval=0.0;
    for(ui d=0;d<dim;d++) ginv(s,d,particles[i].xp)*(particles[i].x[d]-particles[i].xp[d]);
    return retval-2.0*pow(integrator.h)/particles[i].m*particles[i].F[s];
}

template<ui dim> void mpmd<dim>::phet(ui i,ldf eps[dim])
{
    ldf temp[dim];
    for(ui d=0;d<dim;d++) temp[d]=0.0;
    for(ui d=0;d<dim;d++) for(ui s=0;s<dim;s++) for(ui i=0;i<dim;i++) for(ui j=0;j<dim;j++) temp[d]+=ginv(s,d,particles[i].x)*dg(s,i,j,particles[i].x)*eps[i]*eps[j];
    for(ui d=0;d<dim;d++) eps[d]=temp[d]+gimmel(i,d);
}

template<ui dim> void mpmd<dim>::thread_zuiden_wfi(ui i)
{
    ldf eps;
    for(ui d=0;d<dim;d++)
    {
        eps=0.0;
        for(ui dd=0;dd<dim;dd++) eps+=patch.ginv(dd,d,particles[i].xp)*gimmel(i,dd);
        particles[i].xp[d]=particles[i].x[d];
        particles[i].x[d]+=eps;
        particles[i].dx[d]=eps/integrator.h;
    }
}

template<ui dim> void mpmd<dim>::thread_zuiden_protect(ui i)
{
    ui count=0;
    ldf eps[dim],epsp[dim],fit;
    for(ui d=0;d<dim;d++) eps[d]=gimmel(i,d),epsp[d]=eps[d]+1.0;
    do
    {
        count++;
        phet(i,eps);
        fit=0.0;
        for(ui d=0;d<dim;d++) fit+=pow(eps[d]-epsp[d],2),epsp[d]=eps[d];
        if(integrator.generations<=count) break;
    }
    while(fit>numeric_limits<ldf>::epsilon());
    for(ui d=0;d<dim;d++)
    {
        particles[i].xp[d]=particles[i].x[d];
        particles[i].x[d]+=eps[d];
        particles[i].dx[d]=eps[d]/integrator.h;
    }
}

template<ui dim> void mpmd<dim>::thread_zuiden(ui i)
{
    ldf eps[dim],epsp[dim],fit;
    for(ui d=0;d<dim;d++) eps[d]=gimmel(i,d),epsp[d]=eps[d]+1.0;
    do
    {
        phet(i,eps);
        fit=0.0;
        for(ui d=0;d<dim;d++) fit+=pow(eps[d]-epsp[d],2),epsp[d]=eps[d];
    }
    while(fit>numeric_limits<ldf>::epsilon());
    for(ui d=0;d<dim;d++)
    {
        particles[i].xp[d]=particles[i].x[d];
        particles[i].x[d]+=eps[d];
        particles[i].dx[d]=eps[d]/integrator.h;
    }
}

template<ui dim> void mpmd<dim>::thread_calc_forces(ui i)
{
    for(ui j=network.skins[i].size()-1;j<numeric_limits<ui>::max();j--) if(i>network.skins[i][j].neighbor)
    {
        const ldf rsq=embedded_distsq(i,network.skins[i][j].neighbor);
        if(rsq<network.rcosq)
        {
            const ldf r=sqrt(rsq);
            const ldf dVdr=v.dr(network.library[network.skins[i][j].interaction].potential,r,rsq,&network.library[network.skins[i][j].interaction].parameters);
            for(ui d=0;d<dim;d++)
            {
                const ldf delta=embedded_dd(d,i,network.skins[i][j].neighbor);
                const ldf F=delta*dVdr/r;
                #ifdef THREADS
                lock_guard<mutex> freeze(parallel.lock);
                particles[i].F[d]+=F;
                particles[network.skins[i][j].neighbor].F[d]-=F;
                #elif OPENMP
                #pragma omp atomic
                particles[i].F[d]+=F;
                #pragma omp atomic
                particles[network.skins[i][j].neighbor].F[d]-=F;
                #else
                particles[i].F[d]+=F;
                particles[network.skins[i][j].neighbor].F[d]-=F;
                #endif
            }
        }
    }
}

template<ui dim> void mpmd<dim>::integrate()
{
    switch(integrator.method)
    {
        case 2:
            #ifdef THREADS
            for(ui t=0;t<parallel.nothreads;t++) parallel.block[t]=thread([=](ui t){for(ui i=t;i<N;i+=parallel.nothreads) if(!particles[i].fix) thread_zuiden_wfi(i);},t);
            for(ui t=0;t<parallel.nothreads;t++) parallel.block[t].join();
            #elif OPENMP
            #pragma omp parallel for
            for(ui i=0;i<N;i++) if(!particles[i].fix) thread_zuiden_wfi(i);
            #else
            for(ui i=0;i<N;i++) if(!particles[i].fix) thread_zuiden_wfi(i);
            #endif
        break;
        case 1:
            #ifdef THREADS
            for(ui t=0;t<parallel.nothreads;t++) parallel.block[t]=thread([=](ui t){for(ui i=t;i<N;i+=parallel.nothreads) if(!particles[i].fix) thread_zuiden_protect(i);},t);
            for(ui t=0;t<parallel.nothreads;t++) parallel.block[t].join();
            #elif OPENMP
            #pragma omp parallel for
            for(ui i=0;i<N;i++) if(!particles[i].fix) thread_zuiden_protect(i);
            #else
            for(ui i=0;i<N;i++) if(!particles[i].fix) thread_zuiden_protect(i);
            #endif
        break;
        default:
            #ifdef THREADS
            for(ui t=0;t<parallel.nothreads;t++) parallel.block[t]=thread([=](ui t){for(ui i=t;i<N;i+=parallel.nothreads) if(!particles[i].fix) thread_zuiden(i);},t);
            for(ui t=0;t<parallel.nothreads;t++) parallel.block[t].join();
            #elif OPENMP
            #pragma omp parallel for
            for(ui i=0;i<N;i++) if(!particles[i].fix) thread_zuiden(i);
            #else
            for(ui i=0;i<N;i++) if(!particles[i].fix) thread_zuiden(i);
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
