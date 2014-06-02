#ifndef rtgroups_h
#define rtgroups_h

unsigned long long rseedb;

ui randnrb()
{   return rseedb = (16807 * rseedb) % 2147483647;
}

long long rseed;

ldf randnr()
{   return ((rseed = (16807 * rseed) % 2147483647) + .5) / 2147483647.0;
}

bool coinflip()
{   return randnrb() & 32;
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

#endif
