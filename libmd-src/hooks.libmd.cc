#define __libmd_src_file__
#ifndef libmd_h
#include "../libmd.h"
#endif

template<ui dim> ui t_hook<dim>::add(hookptr<dim> p)
{
    //!
    //! This function allows the user to add an userdefined hook which is pointed at by <tt>p</tt>.
    //!
    hooks.push_back(p);
    return hooks.size()-1;
}

template<ui dim> void t_hook<dim>::operator()(ui idx,std::vector<ldf> &parameters,void *sys)
{
    //!
    //! This function calculates a certain hook <tt>hooks[idx]</tt> <br>
    //! The sys pointer is typically a void pointer to the md or mpmd system (which is cast back by using the macro SYS).
    //!
    hooks[idx](parameters,sys);
}

hooktype::hooktype(ui nohook,std::vector<ldf> &param)
{
    //!
    //! This is the hooktype constructor. <br>
    //! It expects the hooktype number: <tt>nohook</tt> which should be aligned with the t_hook<dim>::hooks vector. <br>
    //! Additionally it need parameters for the hook given by <tt>param</tt>. <br>
    //!
    hook=nohook;
    parameters=param;
}
