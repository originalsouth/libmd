#ifndef libmd_h
#include "../libmd.h"
#endif

template<ui dim> indexer<dim>::celldatatype::celldatatype()
{
    nCells = 0;
    IndexDelta = nullptr;
}

template<ui dim> indexer<dim>::celldatatype::~celldatatype()
{
    if (IndexDelta)
        delete[] IndexDelta;
}

template<ui dim> indexer<dim>::indexer()
{
    method=0;
}

/*** Cell algorithm ***/

template<ui dim> void md<dim>::thread_cell (ui i)
{   ui nNeighbors; // Number of neighbors of a cell
    int CellIndices[dim]; // Indices (0 to Q[d]) of cell
    ldf DissqToEdge[dim][3]; // Distance squared from particle to cell edges
    ui d, j, k, particleId, cellId, dissqToCorner, inttype;
    list<ui>::iterator a, b;
    ui NeighboringCells[indexdata.celldata.totNeighbors]; // Cells to check (accounting for boundary conditions)
    ui NeighborIndex[indexdata.celldata.totNeighbors]; // Index (0 to totNeighbors) of neighboring cell
    //vector<list<ui>> Cells(indexdata.celldata.nCells); //Vector for clang++

    // Determine cell indices
    k = i;
    for (d = dim-1; d < numeric_limits<ui>::max(); d--)
    {   CellIndices[d] = k % indexdata.celldata.Q[d];
        k /= indexdata.celldata.Q[d];
    }
    
    // Determine all neighbors
    nNeighbors = 0;
    for (k = 0; k < indexdata.celldata.totNeighbors; k++)
    {   cellId = 0;
        for (d = 0; d < dim && ((simbox.bcond[d] == 1 && indexdata.celldata.Q[d] != 2) || (CellIndices[d]+indexdata.celldata.IndexDelta[k][d] < (int)indexdata.celldata.Q[d] && CellIndices[d]+indexdata.celldata.IndexDelta[k][d] >= 0)); d++)
            cellId = indexdata.celldata.Q[d] * cellId + (indexdata.celldata.Q[d] + CellIndices[d] + indexdata.celldata.IndexDelta[k][d]) % indexdata.celldata.Q[d];
        if (d == dim)
        {   NeighboringCells[nNeighbors] = cellId;
            NeighborIndex[nNeighbors] = k;
            nNeighbors++;
        }
    }
    
    // Loop over all particles in this cell
    for (a = indexdata.celldata.Cells[i].begin(); a != indexdata.celldata.Cells[i].end(); a++)
    {   particleId = *a;
        for (d = 0; d < dim; d++)
        {   DissqToEdge[d][1] = 0;
            if (indexdata.celldata.Q[d] == 2 && simbox.bcond[d] == 1) // Special case: two cells and pbc
                DissqToEdge[d][0] = DissqToEdge[d][2] = pow(indexdata.celldata.CellSize[d]/2 - fabs(indexdata.celldata.CellSize[d]/2 - fmod(simbox.L[d]/2 + particles[particleId].x[d], indexdata.celldata.CellSize[d])), 2);
            else
            {   DissqToEdge[d][0] = pow(fmod(simbox.L[d]/2 + particles[particleId].x[d], indexdata.celldata.CellSize[d]), 2);
                DissqToEdge[d][2] = pow(indexdata.celldata.CellSize[d] - fmod(simbox.L[d]/2 + particles[particleId].x[d], indexdata.celldata.CellSize[d]), 2);
            }
        }

        // Loop over all remaining particles in the same cell
        for (b = next(a); b != indexdata.celldata.Cells[i].end(); b++)
            if (distsq(particleId, *b) < network.sszsq && network.lookup.count(network.hash(particles[particleId].type, particles[*b].type)))
            {   inttype = network.lookup[network.hash(particles[particleId].type, particles[*b].type)];
                network.skins[particleId].push_back(interactionneighbor(*b, inttype));
                network.skins[*b].push_back(interactionneighbor(particleId, inttype));
            }

        // Loop over neighboring cells
        for (k = 0; k < nNeighbors; k++)
        {   
            // Calculate distance (squared) to closest corner
            dissqToCorner = 0;
            for (d = 0; d < dim; d++)
                dissqToCorner += DissqToEdge[d][indexdata.celldata.IndexDelta[NeighborIndex[k]][d]+1];
            // Ignore cell if it is more than network.sszsq away
            if (dissqToCorner > network.sszsq)
                continue;
            // Check all particles in cell
            j = NeighboringCells[k];
            for (b = indexdata.celldata.Cells[j].begin(); b != indexdata.celldata.Cells[j].end(); b++)
                if (distsq(particleId, *b) < network.sszsq && network.lookup.count(network.hash(particles[particleId].type, particles[*b].type)))
                {   inttype = network.lookup[network.hash(particles[particleId].type, particles[*b].type)];
                    network.skins[particleId].push_back(interactionneighbor(*b, inttype));
                    network.skins[*b].push_back(interactionneighbor(particleId, inttype));
                }
        }
    }
}


template<ui dim> void md<dim>::cell()
{
    ui d, i, k, cellId;
    ldf ssz = sqrt(network.sszsq);
    list<ui>::iterator a, b;
    if (!indexdata.celldata.nCells)
    {   ldf nc = 1;
        for (d = 0; d < dim; d++) nc *= indexdata.celldata.Q[d] = (simbox.L[d] < ssz ? 1 : simbox.L[d]/ssz);
        for (; nc > N; nc /= 2)
        {
            k = 0;
            for (d = 1; d < dim; d++)
                if (indexdata.celldata.Q[k] < indexdata.celldata.Q[d])
                    k = d;
            indexdata.celldata.Q[k] /= 2;
        }
        // Compute and check cell sizes
        for (d = 0; d < dim; d++)
        {
            if (indexdata.celldata.Q[d] < 1)
            {   fprintf(stderr, "Error: Q[%d] should be positive, but is %d!\n", d, indexdata.celldata.Q[d]);
                return;
            }
            if ((indexdata.celldata.CellSize[d] = simbox.L[d]/indexdata.celldata.Q[d]) < ssz && indexdata.celldata.Q[d] > 1)
            {   fprintf(stderr, "Error: Q[%d] is too large! (value = %d, max = %d)\n", d, indexdata.celldata.Q[d], max(1, (int)(simbox.L[d]/ssz)));
                return;
            }
        }
        // Compute nCells and totNeighbors
        indexdata.celldata.nCells = 1;
        indexdata.celldata.totNeighbors = 0;
        for (d = 0; d < dim; d++)
            if (indexdata.celldata.Q[d] > 1) // Ignore dimensions with only one cell
            {   indexdata.celldata.nCells *= indexdata.celldata.Q[d];
                indexdata.celldata.totNeighbors = 3*indexdata.celldata.totNeighbors+1;
            }

        indexdata.celldata.Cells.resize(indexdata.celldata.nCells); //Vector for clang++
        // Declare dynamic arrays
        indexdata.celldata.IndexDelta = new int[indexdata.celldata.totNeighbors][dim]; // Relative position of neighboring cell
        // Determine all (potential) neighbors
        // Start with {0,0,...,0,+1}
        if (indexdata.celldata.totNeighbors > 0)
        {   memset(indexdata.celldata.IndexDelta[0], 0, dim*sizeof(ui));
            for (d = dim-1; indexdata.celldata.Q[d] == 1; d--);
            indexdata.celldata.IndexDelta[0][d] = 1;
            for (i = 1; i < indexdata.celldata.totNeighbors; i++)
            {   memcpy(indexdata.celldata.IndexDelta[i], indexdata.celldata.IndexDelta[i-1], dim*sizeof(ui));
                // Set all trailing +1's to -1
                for (d = dim-1; d < dim && (indexdata.celldata.Q[d] == 1 || indexdata.celldata.IndexDelta[i][d] == 1); d--)
                    if (indexdata.celldata.Q[d] > 1)
                        indexdata.celldata.IndexDelta[i][d] = -1;
                // Increase last not-plus-one by one
                if (d < dim)
                    indexdata.celldata.IndexDelta[i][d]++;
            }
        }
    }
    
    // Put the particles in their cells
    for (i = 0; i < indexdata.celldata.nCells; i++)
    	indexdata.celldata.Cells[i].clear();
    for (i = 0; i < N; i++)
    {   cellId = 0;
        for (d = 0; d < dim; d++)
            cellId = indexdata.celldata.Q[d] * cellId + (ui)((simbox.L[d]/2 + particles[i].x[d]) / indexdata.celldata.CellSize[d]);
        indexdata.celldata.Cells[cellId].push_back(i);
    }

    for (i = 0; i < N; i++)
        network.skins[i].clear();

    #ifdef THREADS
    for(ui t=0;t<parallel.nothreads;t++) parallel.block[t]=thread([=](ui t){for(ui i=t;i<indexdata.celldata.nCells;i+=parallel.nothreads) thread_cell(i);},t);
    for(ui t=0;t<parallel.nothreads;t++) parallel.block[t].join();
    #elif OPENMP
    #pragma omp parallel for
    for(ui i=0;i<indexdata.celldata.nCells;i++) thread_cell(i);
    #else
    for(ui i=0;i<indexdata.celldata.nCells;i++) thread_cell(i);
    #endif
}

template<ui dim> void md<dim>::bruteforce()
{
    for(ui i=0;i<N;i++)
    {
        network.skins[i].clear();
        for(ui j=0;j<N;j++) if(i!=j and distsq(i,j)<network.sszsq)
        {
            const pair<ui,ui> it=network.hash(particles[i].type,particles[j].type);
            if(network.lookup.count(it))
            {
                interactionneighbor in(j,network.lookup[it]);
                network.skins[i].push_back(in);
            }
        }
    }
}

