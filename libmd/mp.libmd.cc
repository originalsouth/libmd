template<ui dim> mp<dim>::mp()
{
    setmp();
}

template<ui dim> void mp<dim>::setmp(ui i)
{
    switch(i)
    {
        case 1:
            parameters.assign(2,1);
            fmp=&GAUSSIANBUMP<dim>;
            dfmp=&dGAUSSIANBUMP<dim>;
            ddfmp=&ddGAUSSIANBUMP<dim>;
        break;
        default:
            parameters.assign(1,1);
            fmp=&FLATSPACE<dim>;
            dfmp=&dFLATSPACE<dim>;
            ddfmp=&ddFLATSPACE<dim>;
        break;
    }
}

template<ui dim> void mp<dim>::setmp(fmpptr f,dfmpptr df,ddfmpptr ddf)
{
    fmp=f;
    dfmp=df;
    ddfmp=ddf;
}

template<ui dim> ldf mp<dim>::f(ldf x[dim])
{
    return fmp(x,&parameters);
}

template<ui dim> ldf mp<dim>::g(ui i,ui j,ldf x[dim])
{
    const ldf kdel=kdelta(i,j);
    return kdel+dfmp(i,x,&parameters)*dfmp(j,x,&parameters);
}

template<ui dim> ldf mp<dim>::G(ui s,ui i,ui j,ldf x[dim])
{
    return ddfmp(s,i,x,&parameters)*dfmp(j,x,&parameters)+dfmp(i,x,&parameters)*ddfmp(s,j,x,&parameters);
}
