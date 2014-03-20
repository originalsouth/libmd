#ifndef libmd_h
#include "../libmd.h"
#endif

const ldf mxinv_eps = 1e-12; // matrix inversion cutoff

template<ui dim> box<dim>::box()
{
    memset(bcond,0,dim*sizeof(uc));
    for(ui d=0;d<dim;d++) L[d]=10.0;
    for(ui i=0;i<dim;i++) for(ui j=0;j<dim;j++)
    {
        vshear[i][j]=0.0;
        Lshear[i][j]=0.0;
        LshearInv[i][j]=0.0;
    }
    for(ui d=0;d<dim;d++)
    {
        Lshear[d][d]=L[d];
        LshearInv[d][d]=1.0/L[d];
    }
    boxShear=false;
}

template<ui dim> void box<dim>::shear_boundary(ui i, ui j, ldf velocity)
{
    vshear[i][j]=velocity;
    for(ui d=0;d<dim;d++)
    {
        Lshear[d][d]=L[d];
        LshearInv[d][d] = 1.0/L[d];
    }
    boxShear=true;
    bcond[j]=BCOND::BOXSHEAR;
}

/*** Matrix inverse (determinant)  from Thomas ***/
// Returns the determinant of A
// When the determinant is nonzero, B is the inverse
template<ui dim> ldf det (ldf Ain[dim][dim], ldf B[dim][dim])
{ ui i, j, k;
  int sgn = 1;
  ldf d = 1, t1, t2;
  
  // clone Ain so that it is unchanged
  ldf A[dim][dim];
  std::copy(&Ain[0][0], &Ain[0][0]+dim*dim,&A[0][0]);
  
  // Initialize B to identity matrix
  for (i = 0; i < dim; i++)
    for (j = 0; j < dim; j++)
      B[i][j] = (i==j?1:0);
  for (i = 0; i < dim; i++)
  { // Look for largest pivot
    j = i;
    for (k = i+1; k < dim; k++)
      if (fabs(A[j][i]) < fabs(A[k][i]))
        j = k;
    // No nonzero pivot implies singular matrix
    if (fabs(A[j][i]) < mxinv_eps)
      return 0;
    // Swap rows if necessary
    if (j > i)
    { sgn = -sgn; // Switch sign of determinant
      for (k = i; k < dim; k++)
        swap(A[i][k], A[j][k]);
      for (k = 0; k < dim; k++)
        swap(B[i][k], B[j][k]);
    }
    // Reduce other rows
    t1 = A[i][i];
    for (j = 0; j < dim; j++)
      if (j != i)
      { t2 = A[j][i];
        for (k = i; k < dim; k++)
          A[j][k] = (t1 * A[j][k] - t2 * A[i][k]) / d;
        for (k = 0; k < dim; k++)
          B[j][k] = (t1 * B[j][k] - t2 * B[i][k]) / d;
      }
    d = t1;
  }
  for (i = 0; i < dim; i++)
    for (j = 0; j < dim; j++)
      B[i][j] /= d;
  return sgn * d;
}


template<ui dim> void box<dim>::invert_box()
{   
    ldf d = det(Lshear, LshearInv);
    if (fabs(d) < mxinv_eps) {
        // singular matrix
        ERROR("Singular matrix encountered during box matrix inversion.");
        exit(EXIT_FAILURE);
    }
}
