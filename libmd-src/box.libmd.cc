#define __libmd_src_file__
#ifndef libmd_h
#include "../libmd.h"
#endif

const ldf mxinv_eps=sqrt(std::numeric_limits<ldf>::epsilon()); // matrix inversion cutoff

template<ui dim> box<dim>::box()
{
    //!
    //! Constructer for the system box. Assigns memory to the
    //! arrays that store the box size, boundary conditions, and shear
    //! conditions, all of which depend on the dimension \c dim.
    //!
    memset(bcond,0,dim*sizeof(uc));
    for(ui d=0;d<dim;d++) L[d]=10.0;
    memset(vshear,0,dim*dim*sizeof(ldf));
    memset(Lshear,0,dim*dim*sizeof(ldf));
    memset(LshearInv,0,dim*dim*sizeof(ldf));
    for(ui d=0;d<dim;d++) Lshear[d][d]=L[d],LshearInv[d][d]=1.0/L[d];
    useLshear=false;
    pbcond=false;
}

template<ui dim> void box<dim>::shear_boundary(ui i,ui j,ldf velocity)
{
    //!
    //! Apply a dynamic shear to  the boundary normal to spatial dimension \c j by shear speed
    //! \c velocity in direction \c i. Requires \f$i\neq j\f$.
    //!
    if(i==j)
    {
        ERROR("shear velocity must be perpendicular to boundary");
        exit(EXIT_FAILURE);
    }
    vshear[i][j]=velocity;
    for(ui d=0;d<dim;d++) Lshear[d][d]=L[d],LshearInv[d][d]=1.0/L[d];
    useLshear=true;
    bcond[j]=BCOND::BOXSHEAR;
}

template<ui dim> void box<dim>::skew_boundary(ui i, ui j, ldf displacement)
{
    //!
    //! Apply a static shear by displacing the boundary normal to spatial
    //! dimension \c j in direction \c i by an amount \c displacement.
    //! Requires \f$i\neq j\f$. The boundary conditions in dimension j
    //! must be set to BOXSHEAR (set by default) or HARD; particle dynamics do
    //! not respect the shear if PERIODIC is used.
    //!
    if(i==j)
    {
        ERROR("shear displacement must be perpendicular to boundary");
        exit(EXIT_FAILURE);
    }
    memset(Lshear,0,dim*dim*sizeof(ldf));
    for(ui d=0;d<dim;d++) Lshear[d][d]=L[d];
    Lshear[i][j]=displacement;
    useLshear=true;
    bcond[j]=BCOND::BOXSHEAR;
    invert_box();
}


template<ui dim> ldf det (ldf Ain[dim][dim], ldf B[dim][dim])
{
    //!
    //! Return the determinant of the array \c Ain, and store its
    //! inverse in the array \c B if the determinant is nonzero.
    //!
    using namespace std;
    ui i, j, k;
    int sgn = 1;
    ldf d = 1, t1, t2;

    // clone Ain so that it is unchanged
    ldf A[dim][dim];
    std::copy(&Ain[0][0], &Ain[0][0]+dim*dim,&A[0][0]);

    // Initialize B to identity matrix
    memset(B,0,dim*dim*sizeof(ldf));
    for(i=0;i<dim;i++) B[i][i]=1.0;
    for (i = 0; i < dim; i++)
    { // Look for largest pivot
        j = i;
        for (k = i+1; k < dim; k++)
            if (abs(A[j][i]) < abs(A[k][i]))
                j = k;
        // No nonzero pivot implies singular matrix
        if (abs(A[j][i]) < mxinv_eps)
            return 0;
        // Swap rows if necessary
        if (j > i)
        {
            sgn *= -1; // Switch sign of determinant
            std::swap_ranges(A[i]+i, A[i]+dim, A[j]+i);
            std::swap_ranges(B[i], B[i]+dim, B[j]);
        }
        // Reduce other rows
        t1 = A[i][i]/d;
        for (j = 0; j < dim; j++)
            if (j != i)
            {
                t2 = A[j][i]/d;
                for (k = i; k < dim; k++)
                    A[j][k] = t1 * A[j][k] - t2 * A[i][k];
                for (k = 0; k < dim; k++)
                    B[j][k] = t1 * B[j][k] - t2 * B[i][k];
            }
        d *= t1;
    }
    for (i = 0; i < dim; i++)
        for (j = 0; j < dim; j++)
            B[i][j] /= d;
    return sgn * d;
}


template<ui dim> void box<dim>::invert_box()
{
    using namespace std;
    //!
    //! Invert the system box \c matrix box<dim>::Lshear and store the result
    //! in \c box<dim>::LshearInv.
    //!
    ldf d = det(Lshear, LshearInv);
    if (abs(d) < mxinv_eps) {
        ERROR("singular matrix encountered during box matrix inversion");
        exit(EXIT_FAILURE);
    }
}
