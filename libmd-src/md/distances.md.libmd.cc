#ifndef libmd_h
#include "../../libmd.h"
#endif

template<ui dim> ldf md<dim>::dap(ui d,ldf ad)
{
    //! 
    //! Given a displacement <tt>ad</tt> and corresponding spatial dimension <tt>d</tt>, calculates
    //! the displacement modulo the periodic box size in that dimension.
    //! 
    //! Under periodic boundary conditions, distances are required to lie in the range
    //! <tt>(-simbox.L[d]/2,simbox.L[d]/2)</tt>. If the displacement <tt>dap</tt> is outside
    //! this range, multiples of <tt>simbox.L[d]</tt> are added or subtracted until the displacement
    //! is within these bounds, and the result is returned.
    //!
    ldf da;
    switch(simbox.bcond[d])
    {
        case BCOND::PERIODIC: da=fabs(ad)<0.5*simbox.L[d]?ad:ad-fabs(ad+0.5*simbox.L[d])+fabs(ad-0.5*simbox.L[d]); break;
        default: da=ad; break;
    }
    return da;
}

template<ui dim> ldf md<dim>::distsq(ui p1,ui p2)
{
    //!
    //! Returns the square of the distance between points indexed by
    //! <tt>p1</tt> and <tt>p2</tt>. If periodic boundary conditions
    //! are used, the distance is between the closest periodic images of 
    //! the two points.
    //!
    //!
    ldf retval=0.0;
    for(ui d=0;d<dim;d++) retval+=pow(dd(d,p1,p2),2);
    return retval;
}

template<ui dim> ldf md<dim>::distsq(ldf x1[dim],ldf x2[dim])
{
    //! \overload
    //!
    //! This version accepts as inputs the coordinates of two points
    //! in arrays <tt>x1[dim]</tt> and <tt>x2[dim]</tt> and returns the squared
    //! (periodic) distance between them.
    //! 
    ldf retval=0.0;
    for(ui d=0;d<dim;d++) retval+=pow(dd(d,x1,x2),2);
    return retval;
}

template<ui dim> ldf md<dim>::distsq(ui p1,ldf x2[dim])
{   
    //! \overload
    //!
    //! This version accepts as inputs a point index <tt>p1</tt> and point 
    //! coordinate array <tt>x2[dim]</tt> and returns the squared
    //! (periodic) distance between them.
    //! 
    ldf retval=0.0;
    for(ui d=0;d<dim;d++) retval+=pow(dd(d,p1,x2),2);
    return retval;
}

template<ui dim> ldf md<dim>::distsq(ldf x1[dim],ui p2)
{   
    //! \overload
    //!
    //! This version accepts as inputs a point 
    //! coordinate array <tt>x1[dim]</tt> and a point index <tt>p2</tt> and returns the squared
    //! (periodic) distance between them.
    //! 
    ldf retval=0.0;
    for(ui d=0;d<dim;d++) retval+=pow(dd(d,x1,p2),2);
    return retval;
}

template<ui dim> ldf md<dim>::dd(ui d,ui p1,ui p2)
{
    //!
    //! Returns the distance between points indexed by
    //! <tt>p1</tt> and <tt>p2</tt> along spatial dimension <tt>d</tt<. 
    //! If periodic boundary conditions
    //! are used, the distance is between the closest periodic images of 
    //! the two points.
    //!
    return dd(d,particles[p1].x,particles[p2].x);
}

template<ui dim> ldf md<dim>::dd(ui d,ldf x1[dim],ldf x2[dim])
{   
    //! \overload
    //!
    //! This version accepts as inputs the coordinates of two points
    //! in arrays <tt>x1[dim]</tt> and <tt>x2[dim]</tt> and returns the
    //! (periodic) distance between them along spatial dimension <tt>d</tt>.
    //! 
    ldf ddd=0;
    if (simbox.boxShear) for(ui mu=0;mu<dim;mu++) // use box matrix to calculate distances
    {
       ldf s=0;
       for(ui nu=0;nu<dim;nu++) s+=simbox.LshearInv[mu][nu]*(x2[nu]-x1[nu]);
       if (simbox.bcond[mu]==BCOND::PERIODIC or simbox.bcond[mu]==BCOND::BOXSHEAR) s=fabs(s)<0.5?s:s-fabs(s+0.5)+fabs(s-0.5);
       ddd += simbox.Lshear[d][mu]*s;
    }
    else
    {
        ldf ad=x2[d]-x1[d];
        ddd=dap(d,ad);
    }
    return ddd;
}

template<ui dim> ldf md<dim>::dd(ui d,ui p1,ldf x2[dim])
{   
    //! \overload
    //!
    //! This version accepts as inputs a point index <tt>p1</tt> and point 
    //! coordinate array <tt>x2[dim]</tt> and returns the
    //! (periodic) distance between them along spatial dimension <tt>d</tt>.
    //! 
    return dd(d,particles[p1].x,x2);
}

template<ui dim> ldf md<dim>::dd(ui d,ldf x1[dim],ui p2)
{   
    //! \overload
    //!
    //! This version accepts as inputs a point 
    //! coordinate array <tt>x1[dim]</tt> and a point index <tt>p2</tt> and returns the 
    //! (periodic) distance between them along spatial dimension <tt>d</tt>.
    //! 
    return dd(d,x1,particles[p2].x);
}

template<ui dim> ldf md<dim>::dv(ui d,ui p1,ui p2)
{   
    //! 
    //! Returns the velocity difference between the closest periodic image 
    //! of two particles indexed by <tt>p1</tt> and <tt>p2</tt>. If the vector
    //! connecting the periodic images crosses a sheared boundary, the resulting
    //! velocity difference takes into account the shear velocity of the periodic
    //! image as well. The result is returned as a \f$d\f$-dimensional array
    //! of velocity components.
    //! 
    ldf dv = particles[p2].dx[d]-particles[p1].dx[d];
    if (simbox.boxShear)
    {
        // use box matrix to calculate boundary crossings, and adjust relative velocity accordingly
        for (ui j=0;j<dim;j++)
        {
            ldf s=0,bcross=0;
            for (ui k=0;k<dim;k++)
            {
               s+=simbox.LshearInv[j][k]*(particles[p2].x[k]-particles[p1].x[k]);
            }
            if (fabs(s) > 0.5 and (simbox.bcond[j]==BCOND::PERIODIC or simbox.bcond[j]==BCOND::BOXSHEAR)) bcross = floor(s+.5);
            dv -= simbox.vshear[d][j]*bcross;
        }
    }
    return dv;
}
