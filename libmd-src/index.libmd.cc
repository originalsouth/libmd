#define __libmd_src_file__
#ifndef libmd_h
#include "../libmd.h"
#endif

template<ui dim> indexer<dim>::celldatatype::celldatatype()
{
    //!
    //! Constructor for <tt>celldatatype</tt>, used by the indexing algorithm md<dim>::cell().
    //!
    IndexDelta = nullptr;
    Cells=nullptr;
    OutsideBox=nullptr;
}

template<ui dim> indexer<dim>::celldatatype::~celldatatype()
{
    //!
    //! Destructor for <tt>celldatatype</tt>, used by the indexing algorithm md<dim>::cell().
    //!
    if (IndexDelta)
    {   delete[] IndexDelta;
        delete[] Cells;
        delete[] OutsideBox;
    }
}

template<ui dim> indexer<dim>::kdtreedatatype::kdtreedatatype()
{
    //!
    //! Constructor for <tt>kdtreedatatype</tt>, used by the indexing algorithm md<dim>::kdtree().
    //!
    Idx = nullptr;
    Pmin = nullptr;
    Pmax = nullptr;
}

template<ui dim> indexer<dim>::kdtreedatatype::~kdtreedatatype()
{
    //!
    //! Destructor for <tt>kdtreedatatype</tt>, used by the indexing algorithm md<dim>::kdtree().
    //!
    if (Idx)
    {   delete[] Idx;
        delete[] Pmin;
        delete[] Pmax;
    }
}

template<ui dim> indexer<dim>::indexer()
{
    //!
    //! Constructor for indexer struct. Initializes <tt>method</tt> to <tt>INDEX::CELL</tt> (md<dim>::cell()).
    //!
    method=INDEX::CELL;
}

template<ui dim> ldf dotprod (ldf A[], ldf B[])
{
    //!
    //! Computes the dot product of <tt>A[]</tt> and <tt>B[]</tt>, assumed to both be of size <tt>dim</tt>.
    //!
    ldf x = 0;
    for (ui d = 0; d < dim; d++)
        x += A[d]*B[d];
    return x;
}
