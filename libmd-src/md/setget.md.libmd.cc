#ifndef libmd_h
#include "../../libmd.h"
#endif

template<ui dim> void md<dim>::set_damping(ldf coefficient)
{
    if(avars.noftypedamping==numeric_limits<ui>::max())
    {
        DEBUG_2("activating damping with damping coefficient: %Lf",coefficient);
        vector<ldf> parameters(1,coefficient);
        avars.noftypedamping=add_forcetype(EXTFORCE::DAMPING,nullptr,&parameters);
        assign_all_forcetype(avars.noftypedamping);
        return true;
    }
    else
    {
        DEBUG_2("modifying damping coefficient to: %Lf",coefficient);
        network.forcelibrary[avars.noftypedamping].parameters[0]=coefficient;       
    }
}

template<ui dim> bool md<dim>::unset_damping()
{
    if(avars.noftypedamping<numeric_limits<ui>::max()) 
    {
        DEBUG_2("disabling damping");
        rem_forcetype(avars.noftypedamping);
        return true;
    }
    else
    {
        WARNING("no damping set");
        return false;
    }
}

template<ui dim> void md<dim>::set_rco(ldf rco)
{
    network.rco=rco;
}

template<ui dim> void md<dim>::set_ssz(ldf ssz)
{
    network.ssz=ssz;
    set_reserve(ssz);
}

template<ui dim> void md<dim>::set_reserve(ldf ssz)
{
    set_reserve(ssz,N);
}

template<ui dim> void md<dim>::set_reserve(ldf ssz,ui M)
{
    ldf area=1.0;
    for(ui d=0;d<dim;d++) area*=simbox.L[d];
    const ldf vol=(pow(acos(-1.0),((ldf)dim)/2.0)/tgamma(1.0+((ldf)dim)/2.0))*pow(ssz,dim);
    const ui fed=(ui)(((ldf)M)*(vol)/(area))*2+4;
    DEBUG_1("reserved skin size set to %u skins",fed);
    for(ui i=0;i<N;i++) network.skins[i].reserve(fed);
}

template<ui dim> void md<dim>::set_type(ui p, ui newtype)
{
    particles[p].type=newtype;
}

template<ui dim> void md<dim>::set_index_method(ui method)
{
    indexdata.method=method;
}
