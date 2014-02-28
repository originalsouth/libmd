#ifndef libmd_h
#include "../../libmd.h"
#endif

template<ui dim> ldf md<dim>::dap(ui i,ldf ad)
{
    ldf d;
    switch(simbox.bcond[i])
    {
        case BCOND::PERIODIC: d=fabs(ad)<0.5*simbox.L[i]?ad:ad-fabs(ad+0.5*simbox.L[i])+fabs(ad-0.5*simbox.L[i]); break;
        default: d=ad; break;
    }
    return d;
}

template<ui dim> ldf md<dim>::distsq(ui p1,ui p2)
{
    ldf retval=0.0;
    for(ui i=0;i<dim;i++) retval+=pow(dd(i,p1,p2),2);
    return retval;
}

template<ui dim> ldf md<dim>::dd(ui i,ui p1,ui p2) //TODO: fix non-periodic boundary conditions plus shear
{
    ldf d=0;
    if (simbox.boxShear)
    {
        // use box matrix to calculate distances
        ldf s;
        for (ui j=0;j<dim;j++)
        {
           s=0;
           for (ui k=0;k<dim;k++)
           {
               s+=simbox.LshearInv[j][k]*(particles[p2].x[k]-particles[p1].x[k]);
           }
           if (simbox.bcond[j]==BCOND::PERIODIC or simbox.bcond[j]==BCOND::BOXSHEAR) s=fabs(s)<0.5?s:s-fabs(s+0.5)+fabs(s-0.5);
           d += simbox.Lshear[i][j]*s;
        }
    }
    else
    {
        ldf ad=particles[p2].x[i]-particles[p1].x[i];
        d=dap(i,ad);
    }
    return d;
}
