#ifndef libmd_h
#include "../libmd.h"
#endif

forcetype::forcetype(ui noexternalforce,vector<vector<ui>> *plist,vector<ldf> *param)
{
    //!
    //! This is the forcetype constructor. <br>
    //! It expects the externalforce number: <tt>noexternalforce</tt> which should be aligned with the externalforces<dim>::extforces vector. <br>
    //! Additionally it need parameters for the external force given by <tt>param</tt>. <br>
    //! Optionally an interacting particle list can be given.
    //!
    externalforce=noexternalforce;
    if(plist) particles=*plist;
    else particles.resize(0);
    parameters=*param;
}

template<ui dim> externalforces<dim>::externalforces()
{
    //!
    //! This is the externalforces<dim> constuctor. <br>
    //! It reserves 8 slots in the extforces vector and adds the builtin externalforces. <br>
    //!
    extforces.reserve(8);
    add(DAMPING<dim>);
    add(DISSIPATION<dim>);
}

template<ui dim> ui externalforces<dim>::add(extforceptr<dim> p)
{
    //!
    //! This function allows the user to add an userdefined external force which is pointed at by <tt>p</tt>.
    //!
    extforces.push_back(p);
    return extforces.size()-1;
}

template<ui dim> void externalforces<dim>::operator()(ui type,ui i,vector<ui> *particles,vector<ldf> *parameters,void *sys)
{
    //!
    //! This function calculates a certain external force <tt>extforces[type]</tt> for particle <tt>i</tt> with interacting particle list <tt>particles</tt>. <br>
    //! The sys pointer is typically a void pointer to the md or mpmd system (which is cast back by using the macro SYS).
    //!
    extforces[type](i,particles,parameters,sys);
}
