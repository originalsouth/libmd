#define __libmd_src_file__
#ifndef libmd_h
#include "../../libmd.h"
#endif

template<ui dim> ui md<dim>::add_hook(ui hook,std::vector<ldf> &parameters)
{
    //!
    //! This function adds a new hook, of the given type and with the given parameters, to <tt>hooks.hookers[]</tt> and returns its index.
    //!
    hooktype temp(hook,parameters);
    hooks.hookers.push_back(temp);
    return hooks.hookers.size()-1;
}

template<ui dim> bool md<dim>::mod_hook(ui htype,ui hook,std::vector<ldf> &parameters)
{
    //!
    //! This function replaces the forcetype in <tt>hook.hookers[]</tt> with index <tt>htype</tt>
    //! with a forcetype of the given type and with the given parameters.
    //! It returns whether the given forcelibrary element exists.
    //!
    if(htype<hooks.hookers.size())
    {
        hooktype temp(hook,parameters);
        hooks.hookers[htype]=temp;
        return true;
    }
    else return false;
}

template<ui dim> bool md<dim>::rem_hook(ui htype)
{
    //!
    //! This function removes the forcetype with index <tt>htype</tt> from <tt>hook.hookers[]</tt>.
    //! It returns whether the given forcelibrary element existed.<br>
    //! Note: the last forcetype in <tt>hook.hookers[]</tt> takes the place of the old one.
    //!
    ui pos=hooks.hookers.size();
    if(htype<pos)
    {
        if(htype<pos-1) hooks.hookers[htype]=hooks.hookers.back();
        hooks.hookers.pop_back();
        return true;
    }
    else return false;
}

template<ui dim> bool md<dim>::run_hook(ui htype)
{
    //!
    //! This function runs the hook with index <tt>htype</tt> from <tt>hook.hookers[]</tt>.
    //!
    if(htype<hooks.hookers.size())
    {
        hooks.hook(hooks.hookers[htype].hook,hooks.hookers[htype].parameters,this);
        return true;
    }
    else return false;
}

template<ui dim> void md<dim>::run_hooks()
{
    //!
    //! This function runs all the hook in <tt>hook.hookers[]</tt>.
    //!
    for(hooktype htype:hooks.hookers) hooks.hook(htype.hook,htype.parameters,this);
}
