#ifndef libmd_h
#include "../libmd.h"
#endif

template<ui dim> indexer<dim>::celldatatype::celldatatype()
{
    IndexDelta = nullptr;
}

template<ui dim> indexer<dim>::celldatatype::~celldatatype()
{
    if (IndexDelta)
        delete[] IndexDelta;
}

template<ui dim> indexer<dim>::kdtreedatatype::kdtreedatatype()
{
    Idx = nullptr;
    Pmin = nullptr;
    Pmax = nullptr;
}

template<ui dim> indexer<dim>::kdtreedatatype::~kdtreedatatype()
{
    if (Idx)
    {   delete[] Idx;
        delete[] Pmin;
        delete[] Pmax;
    }
}

template<ui dim> indexer<dim>::indexer()
{
    method=0;
}

template<ui dim> ldf dotprod (ldf A[], ldf B[])
{   ldf x = 0;
    for (ui d = 0; d < dim; d++)
        x += A[d]*B[d];
    return x;
}
