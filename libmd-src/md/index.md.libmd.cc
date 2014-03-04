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
        case INDEX::KD_TREE:
            kdtree();
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
