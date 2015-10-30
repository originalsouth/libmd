#define __libmd_src_file__
#ifndef libmd_h
#include "../libmd.h"
#endif

template<ui dim> bcond<dim>::bcond()
{
    //!
    //! Constructor of the bcond structure
    //!
    //! This functions reserves space for the vector of functionpointers,
    //! and adds the BCOND_... functions in the order defined by the BCOND enum structure.
    //!
    bcond_p.reserve(8);
    bcond_x.reserve(8);
    add(BCOND_NONE<dim>,BCOND_NONE<dim>);
    add(BCOND_PERIODIC<dim>,BCOND_PERIODIC<dim>);
    add(BCOND_HARD<dim>,BCOND_HARD<dim>);
    add(BCOND_BOXSHEAR<dim>,BCOND_BOXSHEAR<dim>);
}

template<ui dim> ui bcond<dim>::add(bcondpptr<dim> p,bcondxptr<dim> x)
{
    //!
    //! This function allows the user to add (costum) boundary condition functions.
    //!
    //! The first argument <tt>bcondpptr<dim> p</tt> is for particles and the second
    //! argument <tt>bcondxptr<dim> x</tt> is for points (or positions).
    //!
    bcond_p.push_back(p);
    bcond_x.push_back(x);
    return bcond_p.size();
}

template<ui dim> void bcond<dim>::operator()(ui d,ui i,void *sys)
{
    //!
    //! This function applies the by the user defined boundary conditions for certain dimension in simbox.bcond[d]
    //! to a certain particle <tt>i</tt> in <tt>md<dim</tt> system <tt>sys</tt>.
    //!
    bcond_p[SYS->simbox.bcond[d]](d,i,sys);
}

template<ui dim> void bcond<dim>::operator()(ui d,ldf x[dim],void *sys)
{
    //!
    //! This function applies the by the user defined boundary conditions for certain dimension in simbox.bcond[d]
    //! to a certain point <tt>x</tt> in <tt>md<dim</tt> system <tt>sys</tt>.
    //!
    bcond_x[SYS->simbox.bcond[d]](d,x,sys);
}

template<ui dim> void bcond<dim>::operator()(ui k,ui d,ui i,void *sys)
{
    //!
    //! This function applies the by the user invoked boundary conditions for certain dimension
    //! to a certain particle <tt>i</tt> in <tt>md<dim</tt> system <tt>sys</tt>.
    //!
    bcond_p[k](d,i,sys);
}

template<ui dim> void bcond<dim>::operator()(ui k,ui d,ldf x[dim],void *sys)
{
    //!
    //! This function applies the by the user invoked boundary conditions for certain dimension
    //! to a certain point <tt>x</tt> in <tt>md<dim</tt> system <tt>sys</tt>.
    //!
    bcond_x[k](d,x,sys);
}
