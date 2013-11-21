///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  Begin of LIBRARY source                                                                                      //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#include "libmd.h"

#include "libmd/potentials.libmd.cc"
#include "libmd/particle.libmd.cc"
#include "libmd/interact.libmd.cc"
#include "libmd/pairpotentials.libmd.cc"
#include "libmd/integrators.libmd.cc"


/*** Cell ***/

ldf square (ldf x)
{	return x*x;
}

template<ui dim> void md<dim>::cell (ui Q[dim])
{	ui nCells; // Total number of cells (= prod(Q))
	ui totNeighbors; // Total number of (potential) neighboring cells to check (= (3^d-1)/2)
	ui nNeighbors; // Number of neighbors of a cell
	int CellIndices[dim]; // Indices (0 to Q[d]) of cell
	ldf CellSize[dim]; // Length of cell in each dimension
	ldf DissqToEdge[dim][3]; // Distance squared from particle to cell edges
	ui d, i, j, k, particleId, cellId, dissqToCorner, inttype;
	list<ui>::iterator a, b;

	// Compute and check cell sizes
	for (d = 0; d < dim; d++)
	{	if (Q[d] < 1)
		{	fprintf(stderr, "Error: Q[%d] should be positive, but is %d!\n", d, Q[d]);
			return;
		}
		if ((CellSize[d] = simbox.L[d]/Q[d]) < rco)
		{	fprintf(stderr, "Error: Q[%d] is too large! (value = %d, max = %d)\n", d, Q[d], (ui)(simbox.L[d]/rco));
			return;
		}
	}

	// Compute nCells and totNeighbors
	nCells = 1;
	totNeighbors = 0;
	for (d = 0; d < dim; d++)
		if (Q[d] > 1) // Ignore dimensions with only one cell
		{	nCells *= Q[d];
			totNeighbors = 3*totNeighbors+1;
		}

	// Declare dynamic arrays
	int (*IndexDelta)[dim] = new int[totNeighbors][dim]; // Relative position of neighboring cell
	ui *NeighboringCells = new ui[totNeighbors]; // Cells to check (accounting for boundary conditions)
	ui *NeighborIndex = new ui[totNeighbors]; // Index (0 to totNeighbors) of neighboring cell
	list<ui> *Cells = new list<ui>[nCells];

	// Determine all (potential) neighbors
	// Start with {0,0,...,0,+1}
	if (totNeighbors > 0)
	{	memset(IndexDelta[0], 0, dim*sizeof(ui));
		for (d = dim-1; Q[d] == 1; d--);
		IndexDelta[0][d] = 1;
		for (i = 1; i < totNeighbors; i++)
		{	memcpy(IndexDelta[i], IndexDelta[i-1], dim*sizeof(ui));
			// Set all trailing +1's to -1
			for (d = dim-1; d < dim && (Q[d] == 1 || IndexDelta[i][d] == 1); d--)
				if (Q[d] > 1)
					IndexDelta[i][d] = -1;
			// Increase last not-plus-one by one
			if (d < dim)
				IndexDelta[i][d]++;
		}
	}

	// Put the particles in their cells
	for (i = 0; i < N; i++)
	{	cellId = 0;
		for (d = 0; d < dim; d++)
			cellId = Q[d] * cellId + (ui)((simbox.L[d]/2 + particles[i].x[d]) / CellSize[d]);
		Cells[cellId].push_back(i);
	}

	memset(CellIndices, 0, sizeof(CellIndices));
	for (i = 0; i < nCells; i++)
	{	
		// Determine all neighbors
		nNeighbors = 0;
		for (k = 0; k < totNeighbors; k++)
		{	cellId = 0;
			for (d = 0; d < dim && ((simbox.bcond[d] == 1 && Q[d] != 2) || (CellIndices[d]+IndexDelta[k][d] < Q[d] && CellIndices[d]+IndexDelta[k][d] >= 0)); d++)
				cellId = Q[d] * cellId + (Q[d] + CellIndices[d] + IndexDelta[k][d]) % Q[d];
			if (d == dim)
			{	NeighboringCells[nNeighbors] = cellId;
				NeighborIndex[nNeighbors] = k;
				nNeighbors++;
			}
		}
		
		// Loop over all particles in this cell
		for (a = Cells[i].begin(); a != Cells[i].end(); a++)
		{	particleId = *a;
			for (d = 0; d < dim; d++)
			{	DissqToEdge[d][1] = 0;
				if (Q[d] == 2 && simbox.bcond[d] == 1) // Special case: two cells and pbc
					DissqToEdge[d][0] = DissqToEdge[d][2] = square(CellSize[d]/2 - fabs(CellSize[d]/2 - fmod(simbox.L[d]/2 + particles[particleId].x[d], CellSize[d])));
				else
				{	DissqToEdge[d][0] = square(fmod(simbox.L[d]/2 + particles[particleId].x[d], CellSize[d]));
					DissqToEdge[d][2] = square(CellSize[d] - fmod(simbox.L[d]/2 + particles[particleId].x[d], CellSize[d]));
				}
			}

			// Loop over all remaining particles in the same cell
			for (b = next(a); b != Cells[i].end(); b++)
				if (distsq(&particles[particleId], &particles[*b]) < rcosq)
				{	inttype = network.lookup[network.hash(particles[particleId].type, particles[*b].type)];
					network.skins[particleId].push_back(interactionneighbor(*b, inttype));
					network.skins[*b].push_back(interactionneighbor(particleId, inttype));
				}

			// Loop over neighboring cells
			for (k = 0; k < nNeighbors; k++)
			{	
				// Calculate distance (squared) to closest corner
				dissqToCorner = 0;
				for (d = 0; d < dim; d++)
					dissqToCorner += DissqToEdge[d][IndexDelta[NeighborIndex[k]][d]+1];
				// Ignore cell if it is more than rcosq away
				if (dissqToCorner > rcosq)
					continue;
				// Check all particles in cell
				j = NeighboringCells[k];
				for (b = Cells[j].begin(); b != Cells[j].end(); b++)
					if (distsq(&particles[particleId], &particles[*b]) < rcosq)
					{	inttype = network.lookup[network.hash(particles[particleId].type, particles[*b].type)];
						network.skins[particleId].push_back(interactionneighbor(*b, inttype));
						network.skins[*b].push_back(interactionneighbor(particleId, inttype));
					}
			}
		}
		// Indices of next cell
		for (d = dim-1; d < dim && CellIndices[d] == Q[d]-1; d--)
			CellIndices[d] = 0;
		if (d < dim)
			CellIndices[d]++;
	}

	// Clean up
	delete[] Cells;
	delete[] NeighboringCells;
	delete[] NeighborIndex;
	delete[] IndexDelta;
}


template<ui dim> void md<dim>::index()
{	ui Q[dim], i, d, b;
	double nCells = 1;
	for (i = 0; i < N; i++)
		network.skins[i].clear();
	for (d = 0; d < dim; d++)
		Q[d] = simbox.L[d]/rco;
	for (; nCells > N; nCells /= 2)
	{	b = 0;
		for (d = 1; d < dim; d++)
			if (Q[b] < Q[d])
				b = d;
		Q[b] /= 2;
	}
	cell(Q);
}


template<ui dim> md<dim>::md()
{
    N=0;
}

template<ui dim> md<dim>::md(ui particlenr)
{
    N=particlenr;
    particles.resize(N);
    network.skins.resize(N);
}

template<ui dim> void md<dim>::add_typeinteraction(ui type1,ui type2,ui potential,vector<ldf> *parameters)
{
    pair<ui,ui> id=network.hash(type1,type2);
    interactiontype itype(potential,parameters,v(potential,network.rco,network.rcosq,parameters));
    if(network.lookup.find(id)==network.lookup.end()) network.library.push_back(itype),network.lookup[id]=network.library.size()-1;
    else network.library[network.lookup[id]]=itype;
}

template<ui dim> void md<dim>::mod_typeinteraction(ui type1,ui type2,ui potential,vector<ldf> *parameters)
{
    pair<ui,ui> id=network.hash(type1,type2);
    interactiontype itype(potential,parameters,v(potential,network.rco,network.rcosq,parameters));
    if(network.lookup.find(id)==network.lookup.end()) network.library.push_back(itype),network.lookup[id]=network.library.size()-1;
    else network.library[network.lookup[id]]=itype;
}

template<ui dim> ldf md<dim>::distsq(ui p1,ui p2)
{
    ldf retval=0.0;
    for(ui i=0;i<dim;i++)
    {
        ldf ad=fabs(particles[p2].x[i]-particles[p1].x[i]),d;
        switch(simbox.bcond[i])
        {
            case 1: d=(ad<simbox.L[i]/2.0?ad:simbox.L[i]-ad); break;
            default: d=ad; break;
        }
        retval+=pow(d,2);
    }
    return retval;
}

template<ui dim> ldf md<dim>::dd(ui i,ui p1,ui p2)
{
    ldf ad=particles[p2].x[i]-particles[p1].x[i],d;
    switch(simbox.bcond[i])
    {
        case 1: d=fmod(3.0*simbox.L[i]/2.0+ad,simbox.L[i])-simbox.L[i]/2.0; break;
        default: d=ad; break;
    }
    return d;
}

template<ui dim> void md<dim>::thread_clear_forces(ui i)
{
    for(ui d=0;d<dim;d++) particles[i].F[d]=0.0;
}

//TODO: Implement atomics
//TODO: What if potential is velocity dependent (damping)?
template<ui dim> void md<dim>::thread_calc_forces(ui i)
{
    for(ui j=network.skins[i].size()-1;j<numeric_limits<ui>::max();j--) if(i>network.skins[i][j].neighbor)
    {
        const ldf rsq=distsq(i,network.skins[i][j].neighbor);
        if(rsq<network.rcosq)
        {
            const ldf r=sqrt(rsq);
            const ldf dVdr=v.dr(network.library[network.skins[i][j].interaction].potential,r,rsq,&network.library[network.skins[i][j].interaction].parameters);
            for(ui d=0;d<dim;d++)
            {
                const ldf delta=dd(d,i,network.skins[i][j].neighbor);
                const ldf F=delta*dVdr/r;
                particles[i].F[d]+=F; //FIXME: atomic
                particles[network.skins[i][j].neighbor].F[d]-=F; //FIXME: atomic
            }
        }
    }
}

//FIXME: Thomas
template<ui dim> void md<dim>::index()
{
    for(ui i=0;i<N;i++)
    {
        network.skins[i].clear();
        for(ui j=0;j<N;j++) if(i!=j and distsq(i,j)<network.sszsq)
        {
            interactionneighbor in(j,network.lookup[network.hash(i,j)]);
            network.skins[i].push_back(in);
        }
    }
}

//TODO: Make parallel launcher
template<ui dim> void md<dim>::calc_forces()
{
    for(ui i=0;i<N;i++) thread_clear_forces(i);
    for(ui i=0;i<N;i++) thread_calc_forces(i);
}

template<ui dim> void md<dim>::thread_integrate(ui i)
{
    switch(integrator.method)
    {
        case 1:

        break;
        default:
        const ldf o=integrator.h/particles[i].m;
        for(ui d=0;d<dim;d++)
        {
            particles[i].dx[d]+=o*particles[i].F[d];
            particles[i].x[d]+=integrator.h*particles[i].dx[d];
        }
        break;
    }

}

//TODO: Make parallel launcher
//TODO: Implement boundary conditions
//TODO: Implement masses
template<ui dim> void md<dim>::integrate()
{
    for(ui i=0;i<N;i++) thread_integrate(i);
    for(integrator.generation=1;integrator.generation<integrator.generations;integrator.generation++)
    {
        calc_forces();
        for(ui i=0;i<N;i++) thread_integrate(i);
    }
}

template<ui dim> void md<dim>::timestep()
{
    if(network.update) index();
    calc_forces();
    integrate();
}

template<ui dim> void md<dim>::timesteps(ui k)
{
    for(ui i=0;i<k;i++) timestep();
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  End of LIBRARY source                                                                                        //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
