#ifndef libmd_h
#include "../../libmd.h"
#endif

template<ui dim> ldf md<dim>::dap(ui d,ldf ad)
{
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
    ldf retval=0.0;
    for(ui d=0;d<dim;d++) retval+=pow(dd(d,p1,p2),2);
    return retval;
}

template<ui dim> ldf md<dim>::distsq(ldf x1[dim],ldf x2[dim])
{
    ldf retval=0.0;
    for(ui d=0;d<dim;d++) retval+=pow(dd(d,x1,x2),2);
    return retval;
}

template<ui dim> ldf md<dim>::distsq(ui p1,ldf x2[dim])
{
    ldf retval=0.0;
    for(ui d=0;d<dim;d++) retval+=pow(dd(d,p1,x2),2);
    return retval;
}

template<ui dim> ldf md<dim>::distsq(ldf x1[dim],ui p2)
{
    ldf retval=0.0;
    for(ui d=0;d<dim;d++) retval+=pow(dd(d,x1,p2),2);
    return retval;
}

template<ui dim> ldf md<dim>::dd(ui d,ui p1,ui p2)
{
    ldf ddd=0;
    if (simbox.boxShear)
    {
        // use box matrix to calculate distances
        for(ui mu=0;mu<dim;mu++)
        {
           ldf s=0;
           for(ui nu=0;nu<dim;nu++) s+=simbox.LshearInv[mu][nu]*(particles[p2].x[nu]-particles[p1].x[nu]);
           if (simbox.bcond[mu]==BCOND::PERIODIC or simbox.bcond[mu]==BCOND::BOXSHEAR) s=fabs(s)<0.5?s:s-fabs(s+0.5)+fabs(s-0.5);
           ddd += simbox.Lshear[d][mu]*s;
        }
    }
    else
    {
        ldf ad=particles[p2].x[d]-particles[p1].x[d];
        ddd=dap(d,ad);
    }
    return ddd;
}

template<ui dim> ldf md<dim>::dd(ui d,ldf x1[dim],ldf x2[dim])
{
    ldf ddd=0;
    if (simbox.boxShear)
    {
        // use box matrix to calculate distances
        for(ui mu=0;mu<dim;mu++)
        {
           ldf s=0;
           for(ui nu=0;nu<dim;nu++) s+=simbox.LshearInv[mu][nu]*(x2[nu]-x1[nu]);
           if (simbox.bcond[mu]==BCOND::PERIODIC or simbox.bcond[mu]==BCOND::BOXSHEAR) s=fabs(s)<0.5?s:s-fabs(s+0.5)+fabs(s-0.5);
           ddd += simbox.Lshear[d][mu]*s;
        }
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
    ldf ddd=0;
    if (simbox.boxShear)
    {
        // use box matrix to calculate distances
        for(ui mu=0;mu<dim;mu++)
        {
           ldf s=0;
           for(ui nu=0;nu<dim;nu++) s+=simbox.LshearInv[mu][nu]*(x2[nu]-particles[p1].x[nu]);
           if (simbox.bcond[mu]==BCOND::PERIODIC or simbox.bcond[mu]==BCOND::BOXSHEAR) s=fabs(s)<0.5?s:s-fabs(s+0.5)+fabs(s-0.5);
           ddd += simbox.Lshear[d][mu]*s;
        }
    }
    else
    {
        ldf ad=x2[d]-particles[p1].x[d];
        ddd=dap(d,ad);
    }
    return ddd;
}

template<ui dim> ldf md<dim>::dd(ui d,ldf x1[dim],ui p2)
{
    ldf ddd=0;
    if (simbox.boxShear)
    {
        // use box matrix to calculate distances
        for(ui mu=0;mu<dim;mu++)
        {
           ldf s=0;
           for(ui nu=0;nu<dim;nu++) s+=simbox.LshearInv[mu][nu]*(particles[p2].x[nu]-x1[nu]);
           if (simbox.bcond[mu]==BCOND::PERIODIC or simbox.bcond[mu]==BCOND::BOXSHEAR) s=fabs(s)<0.5?s:s-fabs(s+0.5)+fabs(s-0.5);
           ddd += simbox.Lshear[d][mu]*s;
        }
    }
    else
    {
        ldf ad=particles[p2].x[d]-x1[d];
        ddd=dap(d,ad);
    }
    return ddd;
}

template<ui dim> ldf md<dim>::dv(ui d,ui p1,ui p2)
{
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
