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

template<ui dim> ldf md<dim>::dd(ui i,ui p1,ui p2)
{
    ldf d=0;
    if (simbox.boxShear)
    {
        // use box matrix to calculate distances
        for (ui j=0;j<dim;j++)
        {
           ldf s=0;
           for (ui k=0;k<dim;k++)
           {
               s+=simbox.LshearInv[j][k]*(particles[p2].x[k]-particles[p1].x[k]);
           }
           if (simbox.bcond[j]==BCOND::PERIODIC or simbox.bcond[j]==BCOND::BOXSHEAR)
               s=fabs(s)<0.5?s:s-fabs(s+0.5)+fabs(s-0.5);
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

template<ui dim> ldf md<dim>::dv(ui i,ui p1,ui p2) 
{
    ldf dv = particles[p2].dx[i]-particles[p1].dx[i];
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
            dv -= simbox.vshear[i][j]*bcross;
        }
    }
    return dv;
}
