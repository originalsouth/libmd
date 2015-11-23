#ifndef rtgroups_h
#include "../../rtgroups.h"
#endif

long long hash_skins (vector<vector<interactionneighbor>>& skins)
{   ui i, j, m, n = skins.size();
    long long h = 0, x;
    for (i = 0; i < n; i++)
    {   m = skins[i].size();
        x = 1;
        for (j = 0; j < m; j++)
            x *= 2 * skins[i][j].neighbor + 37;
        h = x * (h ^ (i+73));
    }
    return h;
}

bool test_indexer (bool shear)
{   rseed = rseed_stdval;
    ui D = 2, n = 3000, d, i, j, b;
    ldf Y[D];
    long long h1, h2, h3;
    md<2> sys(n);
    vector<ldf> V = {1.0};
    sys.add_typeinteraction(0,0,0,V);
    sys.simbox.L[0] = 10.0;
    sys.simbox.L[1] = 100.0;
    ui ssz[] = {1,4,7,12,200};
    uc bc[] = {BCOND::NONE, BCOND::PERIODIC};
    for (ui s : ssz)
    for (b = 0; b < 2; b++)
    {   if (shear)
        {   if (b==0)
            {   sys.simbox.skew_boundary(0, 1, .4*sys.simbox.L[0]);
                sys.simbox.bcond[0] = BCOND::NONE;
            }
            else
                sys.simbox.skew_boundary(1, 0, .4*sys.simbox.L[1]);
        }
        else
            sys.simbox.bcond[0] = sys.simbox.bcond[1] = bc[b];
        sys.set_ssz(s);
        for (i = 0; i < n; i++)
        {   if (sys.simbox.useLshear)
            {   for (d = 0; d < D; d++)
                    Y[d] = (sys.simbox.bcond[d] != BCOND::NONE || irand() % 10 ? urand()-.5 : 2*urand()-1);
                for (d = 0; d < D; d++)
                {   sys.particles[i].x[d] = 0;
                    for (j = 0; j < D; j++)
                        sys.particles[i].x[d] += sys.simbox.Lshear[d][j] * Y[j];
                }
            }
            else
                for (d = 0; d < D; d++)
                    sys.particles[i].x[d] = sys.simbox.L[d] * (sys.simbox.bcond[d] != BCOND::NONE || irand() % 10 ? urand()-.5 : 2*urand()-1);
        }
        sys.indexdata.method = INDEX::CELL;
        sys.index();
        h1 = hash_skins(sys.network.skins);
        sys.indexdata.method = INDEX::BRUTE_FORCE;
        sys.index();
        h2 = hash_skins(sys.network.skins);
        if (h1 != h2)
            return false;
        if (!shear)
        {   sys.indexdata.method = INDEX::KD_TREE;
            sys.index();
            h3 = hash_skins(sys.network.skins);
            if (h2 != h3)
                return false;
        }
    }
    return true;
}

bool test_indexer_noshear()
{   if (test_indexer(false))
        test_success;
    else
        test_fail;
}

bool test_indexer_shear()
{   if (test_indexer(true))
        test_success;
    else
        test_fail;
}
