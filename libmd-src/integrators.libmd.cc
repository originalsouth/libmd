#define __libmd_src_file__
#ifndef libmd_h
#include "../libmd.h"
#endif

integrators::integrators()
{
    //!
    //! Integrator constructor <br>
    //! The default integrator is INTEGRATOR::VVERLET <br>
    //! The default timestep is <tt>h=1e-3</tt> <br>
    //! The default number of generations is <tt>generations=10</tt> and only relevant for mpmd<dim>::thread_zuiden_protect. <br>
    //!
    method=INTEGRATOR::VVERLET;
    generations=10;
    h=1e-3;
}
