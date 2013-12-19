#ifndef libmd_h
#include "../libmd.h"
#endif

template<ui dim> box<dim>::box()
{
    for(ui d=0;d<dim;d++) bcond[d]=0;
    for(ui d=0;d<dim;d++) L[d]=10.0;
    for(ui i=0;i<dim;i++) { 
        for(ui j=0;j<dim;j++) {
            vshear[i][j]=0.0; Lshear[i][j]=0.0;
        }
    }
    for(ui d=0;d<dim;d++)  { Lshear[d][d]=L[d]; LshearInv[d][d] = 1./L[d]; }
    LeesEdwards=false;
}

template<ui dim> void box<dim>::shear_boundary(ui i, ui j, ldf velocity)
{
    vshear[i][j]=velocity;
    for(ui d=0;d<dim;d++)  { Lshear[d][d]=L[d]; LshearInv[d][d] = 1./L[d]; }
    LeesEdwards=true;
}

template<ui dim> void box<dim>::invert_box()
{   
    // hack: only correct for simple shear in one direction
    for(ui i=0;i<dim;i++) { 
        for(ui j=0;j<dim;j++) {
            if (i == j) LshearInv[i][j]=1./Lshear[i][j];
            else LshearInv[i][j]=-Lshear[i][j]/(Lshear[i][i]*Lshear[j][j]);
        }
    }
}
