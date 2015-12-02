#ifndef rtgroups_h
#define rtgroups_h

ui rseed, rseed_stdval = 42;

ui irand()
{   return rseed = (16807ll * rseed) % 2147483647;
}

ldf urand()
{   return (irand() + .5) / 2147483647.;
}

bool coinflip()
{   return irand() & 32;
}

bool skins_consistent (vector<vector<interactionneighbor>>& skins)
{   ui n = skins.size(), i, j;
    ui A[n];
    memset(A, 0, sizeof(A));
    for (i = 0; i < n; i++)
        for (j = skins[i].size()-1; j < UI_MAX; j--)
        {   if (skins[i][j].neighbor == i)
                return false;
            A[i] ^= skins[i][j].neighbor+1;
            A[skins[i][j].neighbor] ^= i+1;
        }
    for (i = 0; i < n; i++)
        if (A[i])
            return false;
    return true;
}

struct compare_particles
{	template<ui dim> bool operator() (particle<dim> p1, particle<dim> p2)
	{	ui d;
		for (d = 0; d < dim-1 && p1.x[d] == p2.x[d]; d++);
		return p1.x[d] < p2.x[d];
	}
};

#endif
