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
    ldf add=pow(patch.f(particles[p2].x)-patch.f(particles[p1].x),2);
    return retval+add;
}

template<ui dim> ldf mpmd<dim>::embedded_dd_p1(ui i,ui p1,ui p2)
{
    ldf ad=particles[p2].x[i]-particles[p1].x[i],d;
    switch(simbox.bcond[i])
    {
        case 1: d=fabs(ad)<0.5*simbox.L[i]?ad:ad-fabs(ad+0.5*simbox.L[i])+fabs(ad-0.5*simbox.L[i]); break;
        default: d=ad; break;
    }
    ad=(patch.f(particles[p2].x)-patch.f(particles[p1].x))*patch.df(i,particles[p1].x);
    return -d-ad;
}

template<ui dim> ldf mpmd<dim>::embedded_dd_p2(ui i,ui p1,ui p2)
{
    ldf ad=particles[p2].x[i]-particles[p1].x[i],d;
    switch(simbox.bcond[i])
    {
        case 1: d=fabs(ad)<0.5*simbox.L[i]?ad:ad-fabs(ad+0.5*simbox.L[i])+fabs(ad-0.5*simbox.L[i]); break;
        default: d=ad; break;
    }
    ad=(patch.f(particles[p2].x)-patch.f(particles[p1].x))*patch.df(i,particles[p2].x);
    return d+ad;
}

template<ui dim> void mpmd<dim>::zuiden_C(ui i,ldf C[dim])
{
    ldf temp[dim]={0.0};
    for(ui d=0;d<dim;d++) C[d]=0.0;
    for(ui d=0;d<dim;d++) for(ui mu=0;mu<dim;mu++) temp[d]+=patch.g(d,mu,particles[i].xp)*(particles[i].x[mu]-particles[i].xp[mu]);
    for(ui d=0;d<dim;d++) temp[d]+=pow(integrator.h,2)/particles[i].m*particles[i].F[d];
    for(ui d=0;d<dim;d++) for(ui sigma=0;sigma<dim;sigma++) C[d]+=patch.ginv(d,sigma,particles[i].x)*temp[sigma];
}

template<ui dim> void mpmd<dim>::zuiden_A(ui i,ldf eps[dim])
{
    ldf temp[dim];
    for(ui d=0;d<dim;d++) temp[d]=eps[d];
    for(ui d=0;d<dim;d++) eps[d]=0.0;
    for(ui d=0;d<dim;d++) for(ui mu=0;mu<dim;mu++) for(ui nu=0;nu<dim;nu++) for(ui sigma=0;sigma<dim;sigma++) eps[d]+=patch.ginv(d,sigma,particles[i].x)*patch.dg(sigma,mu,nu,particles[i].x)*temp[mu]*temp[nu];
}

template<ui dim> void mpmd<dim>::thread_zuiden_wfi(ui i)
{
    ldf eps[dim];
    zuiden_C(i,eps);
    for(ui d=0;d<dim;d++) particles[i].xp[d]=particles[i].x[d];
    for(ui d=0;d<dim;d++) particles[i].x[d]+=eps[d];
    for(ui d=0;d<dim;d++) particles[i].dx[d]=eps[d]/integrator.h;
}

template<ui dim> void mpmd<dim>::thread_zuiden_protect(ui i)
{
    ui count=0;
    ldf val,eps[dim],epsp[dim],C[dim];
    zuiden_C(i,C);
    for(ui d=0;d<dim;d++) epsp[d]=eps[d]=C[d];
    do
    {
        val=0.0;
        zuiden_A(i,eps);
        for(ui d=0;d<dim;d++) eps[d]+=C[d];
        for(ui d=0;d<dim;d++) val+=fabs(eps[d]-epsp[d]);
        for(ui d=0;d<dim;d++) epsp[d]=eps[d];
    }
    while(count<integrator.generations or val>numeric_limits<ldf>::epsilon());
    for(ui d=0;d<dim;d++) particles[i].xp[d]=particles[i].x[d];
    for(ui d=0;d<dim;d++) particles[i].x[d]+=eps[d];
    for(ui d=0;d<dim;d++) particles[i].dx[d]=eps[d]/integrator.h;
}

template<ui dim> void mpmd<dim>::thread_zuiden(ui i)
{
    ldf val,eps[dim],epsp[dim],C[dim];
    zuiden_C(i,C);
    for(ui d=0;d<dim;d++) epsp[d]=eps[d]=C[d];
    do
    {
        val=0.0;
        zuiden_A(i,eps);
        for(ui d=0;d<dim;d++) eps[d]+=C[d];
        for(ui d=0;d<dim;d++) val+=fabs(eps[d]-epsp[d]);
        for(ui d=0;d<dim;d++) epsp[d]=eps[d];
    }
    while(val>numeric_limits<ldf>::epsilon());
    for(ui d=0;d<dim;d++) particles[i].xp[d]=particles[i].x[d];
    for(ui d=0;d<dim;d++) particles[i].x[d]+=eps[d];
    for(ui d=0;d<dim;d++) particles[i].dx[d]=eps[d]/integrator.h;
}

template<ui dim> void mpmd<dim>::thread_calc_forces(ui i)
{
    for(ui j=network.skins[i].size()-1;j<numeric_limits<ui>::max();j--) if(i>network.skins[i][j].neighbor)
    {
        const ldf rsq=embedded_distsq(i,network.skins[i][j].neighbor);
        if(rsq<network.rcosq)
        {
            const ldf r=sqrt(rsq);
            const ldf dVdr=v.dr(network.library[network.skins[i][j].interaction].potential,r,&network.library[network.skins[i][j].interaction].parameters);
            for(ui d=0;d<dim;d++)
            {
                #ifdef THREADS
                lock_guard<mutex> freeze(parallel.lock);
                particles[i].F[d]-=embedded_dd_p1(d,i,network.skins[i][j].neighbor)*dVdr/r;
                particles[network.skins[i][j].neighbor].F[d]-=embedded_dd_p2(d,i,network.skins[i][j].neighbor)*dVdr/r;
                #elif OPENMP
                #pragma omp atomic
                particles[i].F[d]-=embedded_dd_p1(d,i,network.skins[i][j].neighbor)*dVdr/r;
                #pragma omp atomic
                particles[network.skins[i][j].neighbor].F[d]-=embedded_dd_p2(d,i,network.skins[i][j].neighbor)*dVdr/r;
                #else
                particles[i].F[d]-=embedded_dd_p1(d,i,network.skins[i][j].neighbor)*dVdr/r;
                particles[network.skins[i][j].neighbor].F[d]-=embedded_dd_p2(d,i,network.skins[i][j].neighbor)*dVdr/r;
                #endif
            }
        }
    }
}

template<ui dim> void mpmd<dim>::integrate()
{
    switch(integrator.method)
    {
        case 4:
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
        case 3:
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
