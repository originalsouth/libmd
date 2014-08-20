#ifndef libmd_h
#include "../../libmd.h"
#endif

template<ui dim> void md<dim>::thread_index_stick(ui i)
{
    //!
    //! This function copies the current position <tt>x[]</tt> of particle <tt>i</tt> to <tt>xsk[]</tt>, for use by md<dim>::test_index().
    //!
    memcpy(particles[i].xsk,particles[i].x,dim*sizeof(ldf));
}

template<ui dim> void md<dim>::index()
{
    //!
    //! This function finds all the pairs of particles that are within <tt>network.ssz</tt> of each other and have an interaction.
    //! They are then put in each other's skinlist (<tt>network.skins[]</tt>).<br>
    //! It calls the algorithm specificied by <tt>indexdata.method</tt>: md<dim>::bruteforce(), md<dim>::kdtree() or md<dim>::cell().
    //!
    for(ui i=0;i<N;i++) thread_index_stick(i);
    switch(indexdata.method)
    {
        case INDEX::BRUTE_FORCE:
            bruteforce();
        break;
        case INDEX::KD_TREE:
            kdtree();
        break;
        default:
            cell();
        break;
    }
    avars.reindex=false;
}

template<ui dim> bool md<dim>::test_index()
{
    //!
    //! This function determines if re-indexing is necessary, based on the position of each particle at the last time of indexing (<tt>xsk[]</tt>):
    //! If no particle has moved by more than \f$(ssz-rco)/2\f$, then there is no pair of particles
    //! that were more than <tt>network.ssz</tt> apart, but now less than <tt>network.rco</tt> apart, so re-indexing is not necessary.
    //! Otherwise, re-indexing might be necessary and thus should be done.
    //!
    ldf delta=pow(network.ssz-network.rco,2)/4.0;
    for(ui i=0;i<N;i++)
    {
        ldf test=0.0;
        for(ui d=0;d<dim;d++) test+=pow(dap(d,particles[i].xsk[d]-particles[i].x[d]),2);
        if(test>delta) return true;
    }
    return false;
}

/*** k-d tree ***/

// Returns the index of the median
template<ui dim> ui md<dim>::kdtree_build (ui first, ui last, ui level)
{   //!
    //! This function recursively builds the k-d tree, to be used by md<dim>::kdtree().<br>
    //! The particles in the subtree marked by <tt>first</tt> and <tt>last</tt> at recursion depth <tt>level</tt>
    //! are split into two subtrees of (about) equal size.
    //! The particles in one subtree have a lower value of their position coordinate in a given dimension than the particles in the other subtree.
    //! The dimension is specified by <tt>indexdata.kdtreedata.DivideByDim[level]</tt>, which is set by md<dim>::kdtree().
    //! The two subtrees are build recursively.
    //!
    if (last - first == 1) // Leaf
    {   // Set minimum and maximum value equal to the position of this particle (Idx[first])
        memcpy(indexdata.kdtreedata.Pmin[first], particles[indexdata.kdtreedata.Idx[first]].x, dim*sizeof(ldf));
        memcpy(indexdata.kdtreedata.Pmax[first], particles[indexdata.kdtreedata.Idx[first]].x, dim*sizeof(ldf));
        return first;
    }
    ui m = (first+last)/2; // Median
    ui sortedDim = indexdata.kdtreedata.DivideByDim[level];
    // Put the index of particle with the median value of coordinate dim in its right place,
    // put all particles with lower values below, and all with higher values above
    nth_element(indexdata.kdtreedata.Idx+first, indexdata.kdtreedata.Idx+m, indexdata.kdtreedata.Idx+last,
                [=](ui i, ui j) -> bool { return particles[i].x[sortedDim]<particles[j].x[sortedDim];});
    // Recursively build subtrees
    ui m1 = kdtree_build(first, m, level+1), m2 = kdtree_build(m, last, level+1), d;
    // Determine minimum and maximum value of each coordinate
    for (d = 0; d < dim; d++)
    {   indexdata.kdtreedata.Pmin[m][d] = min(indexdata.kdtreedata.Pmin[m1][d], indexdata.kdtreedata.Pmin[m2][d]);
        indexdata.kdtreedata.Pmax[m][d] = max(indexdata.kdtreedata.Pmax[m1][d], indexdata.kdtreedata.Pmax[m2][d]);
    }
    return m;
}

template<ui dim> void md<dim>::kdtree_index (ui first1, ui last1, ui first2, ui last2)
{   //!
    //! This function recursively finds neighboring particles, one from subtree 1, the other from subtree 2 (possibly the same subtree).<br>
    //! First it is checked if the subtrees are too far apart; nothing is done then.
    //! Otherwise, the function is called for every combination of subtrees.
    //! When a leaf is reached, the single particle is checked against the particle<b></b>(s) in the other subtree.
    //!
    ui m1 = (first1+last1)/2, m2 = (first2+last2)/2, d;
    ldf sszsq = pow(network.ssz,2);
    // Base cases
    if (m1 == first1 || m2 == first2) // A single particle
    {   // Note: the other subtree contains either one or two particles
        if (m1 != m2 && distsq(indexdata.kdtreedata.Idx[m1], indexdata.kdtreedata.Idx[m2]) < sszsq)
            skinner(indexdata.kdtreedata.Idx[m1], indexdata.kdtreedata.Idx[m2]);
        if (m2 != first2 && distsq(indexdata.kdtreedata.Idx[m1], indexdata.kdtreedata.Idx[first2]) < sszsq)
            skinner(indexdata.kdtreedata.Idx[m1], indexdata.kdtreedata.Idx[first2]);
        if (m1 != first1 && distsq(indexdata.kdtreedata.Idx[first1], indexdata.kdtreedata.Idx[m2]) < sszsq)
            skinner(indexdata.kdtreedata.Idx[first1], indexdata.kdtreedata.Idx[m2]);
        return;
    }
    // Compute distance (squared) between subtrees
    if (m1 != m2) // Note: m1 == m2 iff the subtrees are the same
    {   ldf dissqBetweenSubtrees = 0;
        for (d = 0; d < dim; d++)
        {   if (simbox.bcond[d] == BCOND::PERIODIC)
            {   if (indexdata.kdtreedata.Pmin[m1][d] > indexdata.kdtreedata.Pmax[m2][d])
                    dissqBetweenSubtrees += pow(min(indexdata.kdtreedata.Pmin[m1][d] - indexdata.kdtreedata.Pmax[m2][d],
                                                    simbox.L[d] + indexdata.kdtreedata.Pmin[m2][d] - indexdata.kdtreedata.Pmax[m1][d]), 2);
                else if (indexdata.kdtreedata.Pmin[m2][d] > indexdata.kdtreedata.Pmax[m1][d])
                    dissqBetweenSubtrees += pow(min(indexdata.kdtreedata.Pmin[m2][d] - indexdata.kdtreedata.Pmax[m1][d],
                                                    simbox.L[d] + indexdata.kdtreedata.Pmin[m1][d] - indexdata.kdtreedata.Pmax[m2][d]), 2);
            }
            else if (indexdata.kdtreedata.Pmin[m1][d] > indexdata.kdtreedata.Pmax[m2][d])
                dissqBetweenSubtrees += pow(indexdata.kdtreedata.Pmin[m1][d] - indexdata.kdtreedata.Pmax[m2][d], 2);
            else if (indexdata.kdtreedata.Pmin[m2][d] > indexdata.kdtreedata.Pmax[m1][d])
                dissqBetweenSubtrees += pow(indexdata.kdtreedata.Pmin[m2][d] - indexdata.kdtreedata.Pmax[m1][d], 2);
        }
        if (dissqBetweenSubtrees >= sszsq) // Return if the subtrees are too far apart
            return;
    }
    // Recursively check subtrees
    kdtree_index(first1, m1, first2, m2);
    kdtree_index(first1, m1, m2, last2);
    if (m1 != m2)
        kdtree_index(m1, last1, first2, m2);
    kdtree_index(m1, last1, m2, last2);
}

template<ui dim> void md<dim>::kdtree()
{   //!
    //! This indexing function uses a k-d tree (a multidimensional binary tree).<br>
    //! It makes use of the recursive functions md<dim>::kdtree_build() and md<dim>::kdtree_index().<br>
    //! It does not work if there is both shear and a dimension with periodic boundary conditions.
    //!
    if (simbox.boxShear)
        for (ui d = 0; d < dim; d++)
            if (simbox.bcond[d] == BCOND::PERIODIC)
            {   ERROR("the kd-tree algorithm does not work with both shear and periodic boundary conditions");
                return;
            }
    if (indexdata.kdtreedata.Idx == nullptr || sizeof(indexdata.kdtreedata.Idx) != N*sizeof(ui))
    {   delete[] indexdata.kdtreedata.Idx;
        delete[] indexdata.kdtreedata.Pmin;
        delete[] indexdata.kdtreedata.Pmax;
        indexdata.kdtreedata.Idx = new ui[N];
        indexdata.kdtreedata.Pmin = new ldf[N][dim];
        indexdata.kdtreedata.Pmax = new ldf[N][dim];
    }
    ldf S[dim];
    ui i, n, d, b;
    for (i = 0; i < N; i++)
        indexdata.kdtreedata.Idx[i] = i;
    // Decide on which dimensions to divide the particles by at each recursion level
    memcpy(S, simbox.L, sizeof(S));
    for (i = 0, n = N; n > 1; i++, n = (n+1)/2)
    {   // Look for largest dimension
        b = 0;
        for (d = 1; d < dim; d++)
            if (S[b] < S[d])
                b = d;
        indexdata.kdtreedata.DivideByDim[i] = b;
        S[b] /= 2; // Assume that the system is nicely split into two parts
    }
    kdtree_build(0, N, 0);
    for (i = 0; i < N; i++)
        network.skins[i].clear();
    kdtree_index(0, N, 0, N);
}

/*** Cell algorithm ***/

template<ui dim> void md<dim>::thread_cell (ui c)
{   
    //!
    //! This function finds the neighbors of all particles in cell <tt>c</tt>.
    //! It does this by considering (half of) all neighboring cells, such that no pairs of cells are processed twice.
    //!
    ui nNeighbors; // Number of neighbors of a cell
    int CellIndices[dim]; // Indices (0 to Q[d]) of cell
    ldf DissqToEdge[dim][3]; // Distance squared from particle to cell edges
    ui NeighboringCells[indexdata.celldata.totNeighbors]; // Cells to check (accounting for boundary conditions)
    ui NeighborIndex[indexdata.celldata.totNeighbors]; // Index (0 to totNeighbors) of neighboring cell
    ui d, i, j, k, p1, p2, cellId;
    int ci;
    ldf dissqToCorner, sszsq = pow(network.ssz,2);

    // Determine cell indices
    k = c;
    for (d = dim-1; d < UI_MAX; d--)
    {   DEBUG_3("indexdata.celldata.Q[%u]= %u", d, indexdata.celldata.Q[d]);
        CellIndices[d] = k % indexdata.celldata.Q[d];
        k /= indexdata.celldata.Q[d];
    }

    // Determine all neighbors
    nNeighbors = 0;
    for (k = 0; k < indexdata.celldata.totNeighbors; k++)
    {   cellId = 0;
        for (d = 0; d < dim && (((ci = CellIndices[d]+indexdata.celldata.IndexDelta[k][d]) < (int)indexdata.celldata.Q[d] && ci >= 0)
                                || (simbox.bcond[d] == BCOND::PERIODIC && indexdata.celldata.Q[d] != 2)); d++)
            cellId = indexdata.celldata.Q[d] * cellId + (indexdata.celldata.Q[d] + ci) % indexdata.celldata.Q[d];
        if (d == dim)
        {   NeighboringCells[nNeighbors] = cellId;
            NeighborIndex[nNeighbors] = k;
            nNeighbors++;
        }
    }

    // Loop over all particles in this cell
    for (i = indexdata.celldata.Cells[c].size()-1; i < UI_MAX; i--)
    {   p1 = indexdata.celldata.Cells[c][i];
        for (d = 0; d < dim; d++)
        {   DissqToEdge[d][1] = 0;
            if (indexdata.celldata.Q[d] == 2 && simbox.bcond[d] == BCOND::PERIODIC) // Special case: two cells and pbc
                DissqToEdge[d][0] = DissqToEdge[d][2] = pow(indexdata.celldata.CellSize[d]/2 - fabs((CellIndices[d]+.5) * indexdata.celldata.CellSize[d] - simbox.L[d]/2 - particles[p1].x[d]), 2);
            else
            {   DissqToEdge[d][0] = pow(simbox.L[d]/2 + particles[p1].x[d] - CellIndices[d] * indexdata.celldata.CellSize[d], 2);
                DissqToEdge[d][2] = pow((CellIndices[d]+1) * indexdata.celldata.CellSize[d] - simbox.L[d]/2 - particles[p1].x[d], 2);
            }
        }

        // Loop over all remaining particles in the same cell
        for (j = i-1; j < UI_MAX; j--)
            if (distsq(p1, p2 = indexdata.celldata.Cells[c][j]) < sszsq)
                skinner(p1,p2);

        // Loop over neighboring cells
        for (k = 0; k < nNeighbors; k++)
        {
            // Calculate distance (squared) to closest corner
            dissqToCorner = 0;
            for (d = 0; d < dim; d++)
                dissqToCorner += DissqToEdge[d][indexdata.celldata.IndexDelta[NeighborIndex[k]][d]+1];
            // Ignore cell if it is more than sszsq away
            if (!simbox.boxShear && dissqToCorner > sszsq)
                continue;
            // Check all particles in cell
            for (ui p2 : indexdata.celldata.Cells[NeighboringCells[k]])
                if (distsq(p1,p2) < sszsq)
                    skinner(p1,p2);
        }
    }
}


template<ui dim> void md<dim>::cell()
{   //!
    //! This indexing algorithm uses the cell method.<br>
    //! The system is divided into cells of equal size, which are at least <tt>network.ssz</tt> in size in each dimension.
    //! The particles are put in their corresponding cells.
    //! Two particles that are in cells that do not share at least a corner will be more than <tt>network.ssz</tt> apart,
    //! therefore such pairs of cells need not be checked.<br>
    //! This function does not work if not all the particles are inside the system (as defined by <tt>simbox</tt>).
    //!
    DEBUG_2("exec is here");
    if (network.ssz <= 0)
    {   ERROR("skinsize is not positive (network.ssz = %Lf)", network.ssz);
        return;
    }
    ui c, d, i, k, cellId;
    ldf x;
    list<ui>::iterator a, b;
    ldf nc = 1;
    if (simbox.boxShear)
    {   ldf R;
        for (d = 0; d < dim; d++)
        {   R = pow(dotprod<dim>(simbox.LshearInv[d], simbox.LshearInv[d]), -.5);
            nc *= indexdata.celldata.Q[d] = (R < network.ssz ? 1 : R/network.ssz);
        }
    }
    else
        for (d = 0; d < dim; d++)
            nc *= indexdata.celldata.Q[d] = (simbox.L[d] < network.ssz ? 1 : simbox.L[d]/network.ssz);
    // If number of cells is very large (ssz very small): reduce until number of cells is in the order of N
    for (; nc > N; nc /= 2)
    {
        k = 0;
        for (d = 1; d < dim; d++)
            if (indexdata.celldata.Q[k] < indexdata.celldata.Q[d])
                k = d;
        indexdata.celldata.Q[k] = (indexdata.celldata.Q[k]+1)/2;
    }
    for (d = 0; d < dim; d++)
        DEBUG_3("indexdata.celldata.Q[%u] = %u (from %Lf / %Lf originally)", d, indexdata.celldata.Q[d], simbox.L[d], network.ssz);
    // Compute and check cell sizes
    for (d = 0; d < dim; d++)
        indexdata.celldata.CellSize[d] = simbox.L[d]/indexdata.celldata.Q[d];

    // Compute nCells and totNeighbors
    indexdata.celldata.nCells = 1;
    indexdata.celldata.totNeighbors = 0;
    for (d = 0; d < dim; d++)
        if (indexdata.celldata.Q[d] > 1) // Ignore dimensions with only one cell
        {   indexdata.celldata.nCells *= indexdata.celldata.Q[d];
            indexdata.celldata.totNeighbors = 3*indexdata.celldata.totNeighbors+1;
        }

    // Declare dynamic arrays
    if (indexdata.celldata.IndexDelta == nullptr || sizeof(indexdata.celldata.Cells) != indexdata.celldata.nCells*sizeof(ui)
        || sizeof(indexdata.celldata.IndexDelta) != indexdata.celldata.totNeighbors*dim*sizeof(int))
    {   delete[] indexdata.celldata.Cells;
        delete[] indexdata.celldata.IndexDelta;
        indexdata.celldata.Cells = new vector<ui>[indexdata.celldata.nCells];
        indexdata.celldata.IndexDelta = new int[indexdata.celldata.totNeighbors][dim]; // Relative position of neighboring cell
    }
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

    // Put the particles in their cells
    for (c = 0; c < indexdata.celldata.nCells; c++)
        indexdata.celldata.Cells[c].clear();
    for (i = 0; i < N; i++)
    {   cellId = 0;
        for (d = 0; d < dim; d++)
        {   x = (simbox.boxShear ? dotprod<dim>(simbox.LshearInv[d], particles[i].x) : particles[i].x[d] / simbox.L[d]);
            if (fabs(x) > .5+1e-9)
            {   ERROR("particle %u is outside the simbox: the cell algorithm cannot cope with that",i);
                return;
            }
            cellId = indexdata.celldata.Q[d] * cellId + (x < -.5+3e-9 ? 0 : (ui)(indexdata.celldata.Q[d]*(x+.5-2e-9)));
        }
        indexdata.celldata.Cells[cellId].push_back(i);
    }

    for (i = 0; i < N; i++)
        network.skins[i].clear();

    for(ui c=0;c<indexdata.celldata.nCells;c++) thread_cell(c);
}

template<ui dim> void md<dim>::bruteforce()
{
    //!
    //! This brute force indexing algorithm simply checks every pair of particles.<br>
    //! This function is to be preferred over md<dim>::cell() and md<dim>::kdtree()
    //! if all (or most) particles are within <tt>network.ssz</tt> of each other.
    //!
    DEBUG_2("exec is here");
    for(ui i=0;i<N;i++) network.skins[i].clear();
    ldf sszsq=pow(network.ssz,2);
    for(ui i=0;i<N;i++) for(ui j=i+1;j<N;j++) if(distsq(i,j)<sszsq) skinner(i,j);
}

template<ui dim> void md<dim>::skinner(ui i, ui j)
{
    //!
    //! This function is used by the indexing algorithms.
    //! It checks if particles <tt>i</tt> and <tt>j</tt> (which are within <tt>network.ssz</tt> of each other)
    //! have a (superparticle) type interaction; if so, it adds them (with the interaction) to each other's skinlist.<br>
    //!
    const ui K=network.spid[i];
    pair<ui,ui> it;
    if(K<UI_MAX and K==network.spid[j] and network.sptypes[network.superparticles[K].sptype].splookup.count(it=network.hash(network.superparticles[K].particles[i],network.superparticles[K].particles[j])))
    {
        interactionneighbor in(j,network.sptypes[network.superparticles[K].sptype].splookup[it]);
        network.skins[i].push_back(in);
        in.neighbor=i;
        network.skins[j].push_back(in);
        DEBUG_3("super particle skinned (i,j)=(%u,%u) in %u with interaction %u",i,j,K,network.sptypes[network.superparticles[K].sptype].splookup[it]);
    }
    else
    {
        it=network.hash(particles[i].type,particles[j].type);
        if(network.lookup.count(it))
        {
            interactionneighbor in(j,network.lookup[it]);
            network.skins[i].push_back(in);
            in.neighbor=i;
            network.skins[j].push_back(in);
            DEBUG_3("normally skinned (i,j)=(%u,%u) with interaction %u",i,j,network.lookup[it]);
        }
    }
}
