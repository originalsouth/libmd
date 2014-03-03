long long rseed = 42;

ldf randnr()
{	return ((rseed = (16807 * rseed) % 2147483647) + .5) / 2147483647.0;
}

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
{   ui D = 2, n = 3000, d, i, j;
    ldf Y[D];
    long long h1, h2, h3;
    md<2> sys(n);
    vector<ldf> V = {1.0};
    sys.add_typeinteraction(0,0,0,&V);
    sys.simbox.L[0] = 10.0;
    sys.simbox.L[1] = 100.0;
    ui ssz[] = {1,4,7,12};
    uc bc[] = {BCOND::NONE, BCOND::PERIODIC};
    if (shear)
    {   sys.simbox.boxShear = true;
        sys.simbox.Lshear[0][0] = sys.simbox.L[0];
        sys.simbox.Lshear[1][1] = sys.simbox.L[1];
        sys.simbox.Lshear[0][1] = .4*sys.simbox.L[0];
        sys.simbox.Lshear[1][0] = .4*sys.simbox.L[1];
        sys.simbox.invert_box();
    }
    for (ui& s : ssz)
    for (uc& b : bc)
    {   sys.set_ssz(s);
        sys.simbox.bcond[0] = b;
        sys.simbox.bcond[1] = b;
        for (i = 0; i < n; i++)
        {   if (sys.simbox.boxShear)
            {   for (d = 0; d < D; d++)
                    Y[d] = randnr();
                for (d = 0; d < D; d++)
                {   sys.particles[i].x[d] = 0;
                    for (j = 0; j < D; j++)
                        sys.particles[i].x[d] += sys.simbox.Lshear[d][j] * Y[j];
                }
            }
            else
                for (d = 0; d < D; d++)
                    sys.particles[i].x[d] = sys.simbox.L[d] * randnr();
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
{   printf("%s: %s: ",__FILE__,__FUNCTION__);
    return test_indexer(false);
}

bool test_indexer_shear()
{   printf("%s: %s: ",__FILE__,__FUNCTION__);
    return test_indexer(true);
}
