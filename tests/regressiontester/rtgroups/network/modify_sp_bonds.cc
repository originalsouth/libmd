#ifndef rtgroups_h
#include "../../rtgroups.h"
#endif

bool test_modify_sp_bonds()
{	rseed = rseedb = 42;
	ui runs = 100, n = 100, S = 10, T = 3, nst = 10, nTypes = 30, actions = 10*n;
	ui run, action, d, s, t, m, i, j, a, b;
	ui bruteforce_lookup[n][n];
	ui bruteforce_type_lookup[nTypes][nTypes];
	ui bruteforce_sp_lookup[n][n];
	ui bruteforce_sptype_lookup[T][nst][nst];
	bool sp_pos_used[S][nst];
	bool seen[n];
	pair<ui,ui> id;
	md<2> sys(n);
	sys.simbox.L[0] = sys.simbox.L[1] = 10.0;
	sys.simbox.bcond[0] = sys.simbox.bcond[1] = BCOND::PERIODIC;
	sys.set_rco(4.0);
	sys.set_ssz(5.0);
	vector<ldf> V(1);
	for (run = 0; run < runs; run++)
	{	sys.clear();
		sys.init(n);
		memset(bruteforce_type_lookup, -1, sizeof(bruteforce_type_lookup));
		memset(bruteforce_sptype_lookup, -1, sizeof(bruteforce_sptype_lookup));
		memset(bruteforce_sp_lookup, -1, sizeof(bruteforce_sp_lookup));
		memset(sp_pos_used, false, sizeof(sp_pos_used));
		sys.network.superparticles.resize(S);
		sys.network.sptypes.resize(T);

		// Initialize
		for (i = 0; i < n; i++)
		{	sys.set_type(i, randnrb() % nTypes);
			if (coinflip())
			{ do
				{	s = randnrb() % S;
					j = randnrb() % nst;
				}
				while (sp_pos_used[s][j]);
				sp_pos_used[s][j] = true;
				sys.network.superparticles[s].particles[i] = j;
				sys.network.spid[i] = s;
			}
			for (d = 0; d < 2; d++)
			{	sys.particles[i].x[d] = 10*randnr()-5.0;
				sys.particles[i].dx[d] = randnr()-0.5;
			}
		}
		m = 1;
		for (t = 0; t < T; t++)
		{ sys.add_sptype();
			for (i = 0; i < nst; i++)
				for (j = i+1; j < nst; j++)
				{	V[0] = m;
					sys.add_sp_interaction(t,i,j,0,&V);
					bruteforce_sptype_lookup[t][i][j] = bruteforce_sptype_lookup[t][j][i] = m++;
				}
		}
		for (i = 0; i < nTypes; i++)
			for (j = i; j < nTypes; j++)
				if (coinflip())
				{	V[0] = m;
					sys.add_typeinteraction(i,j,0,&V);
					bruteforce_type_lookup[i][j] = bruteforce_type_lookup[j][i] = m++;
				}
		for (i = 0; i < n; i++)
			for (j = i+1; j < n; j++)
			{	if ((s = sys.network.spid[i]) < (ui)-1 && s == sys.network.spid[j])
					bruteforce_sp_lookup[i][j] = bruteforce_sptype_lookup[sys.network.superparticles[s].sptype][sys.network.superparticles[s].particles[i]][sys.network.superparticles[s].particles[j]];
				bruteforce_lookup[i][j] = bruteforce_type_lookup[sys.particles[i].type][sys.particles[j].type];
			}
		sys.index();

		// Mess around
		for (action = 0; action < actions; action++)
		{	i = randnrb() % n;
			do
				j = randnrb() % n;
			while (i==j);
			if (i>j)
				swap(i,j);
			if ((s = sys.network.spid[i]) < (ui)-1 && s == sys.network.spid[j] && randnrb() % 8 < 7)
			{	id = sys.network.hash(sys.network.superparticles[s].particles[i], sys.network.superparticles[s].particles[j]);
				if (!sys.network.sptypes[sys.network.superparticles[s].sptype].splookup.count(id))
				{	V[0] = m;
					sys.add_sp_bond(i,j,0,&V);
					bruteforce_sp_lookup[i][j] = m++;
				}
				else if (coinflip())
				{	V[0] = m;
					sys.mod_sp_bond(i,j,0,&V);
					bruteforce_sp_lookup[i][j] = m++;
				}
				else
				{	sys.rem_sp_bond(i,j);
					bruteforce_sp_lookup[i][j] = -1;
				}
			}
			else
			{	id = sys.network.hash(sys.particles[i].type, sys.particles[j].type);
				if (!sys.network.lookup.count(id))
				{	V[0] = m;
					sys.add_bond(i,j,0,&V);
					bruteforce_lookup[i][j] = m++;
				}
				else if (coinflip())
				{ V[0] = m;
					sys.mod_bond(i,j,0,&V);
					bruteforce_lookup[i][j] = m++;
				}
				else
				{	sys.rem_bond(i,j);
					bruteforce_lookup[i][j] = -1;
				}
			}
		}

		// Check skin consistency
		if (!skins_consistent(sys.network.skins))
		{	printf("run %d: inconsistency in skins\n", run);
			test_fail;
		}
		// Check lookups
		for (i = 0; i < n; i++)
			for (j = i+1; j < n; j++)
			{	id = sys.network.hash(sys.particles[i].type, sys.particles[j].type);
				if (sys.network.lookup.count(id))
				{	if (sys.network.library[sys.network.lookup[id]].parameters[0] != (ldf)bruteforce_lookup[i][j])
					{	printf("run %d, pair (%d,%d): %d != %d\n", run, i, j, (int)(.5+sys.network.library[sys.network.lookup[id]].parameters[0]), bruteforce_lookup[i][j]);
						test_fail;
					}
				}
				else if (bruteforce_lookup[i][j] != (ui)-1)
				{	printf("run %d, pair (%d,%d): -1 != %d\n", run, i, j, bruteforce_lookup[i][j]);
					test_fail;
				}
				if ((s = sys.network.spid[i]) < (ui)-1 && s == sys.network.spid[j] && sys.network.sptypes[t = sys.network.superparticles[s].sptype].splookup.count(id = sys.network.hash(sys.network.superparticles[s].particles[i], sys.network.superparticles[s].particles[j])))
				{	if (sys.network.library[sys.network.sptypes[t].splookup[id]].parameters[0] != (ldf)bruteforce_sp_lookup[i][j])
					{	printf("run %d, pair (%d,%d) [sp]: %d != %d\n", run, i, j, (int)(.5+sys.network.library[sys.network.sptypes[t].splookup[id]].parameters[0]), bruteforce_sp_lookup[i][j]);
						test_fail;
					}
				}
				else if (bruteforce_sp_lookup[i][j] != (ui)-1)
				{	printf("run %d, pair (%d,%d) [sp]: -1 != %d\n", run, i, j, bruteforce_sp_lookup[i][j]);
					test_fail;
				}
			}
		// Check skins
		for (i = 0; i < n; i++)
		{	memset(seen, false, sizeof(seen));
			for (auto sij : sys.network.skins[i])
			{	j = sij.neighbor;
				a = min(i,j);
				b = i^j^a;
				seen[j] = true;
				m = (bruteforce_sp_lookup[a][b] != (ui)-1 ? bruteforce_sp_lookup[a][b] : bruteforce_lookup[a][b]);
				if (sys.network.library[sij.interaction].parameters[0] != (ldf)(m))
				{	printf("run %d, pair (%d,%d) [skin]: %d != %d\n", run, a, b, (int)(.5+sys.network.library[sij.interaction].parameters[0]), m);
					test_fail;
				}
			}
			ldf sszsq = pow(sys.network.ssz,2);
			for (j = 0; j < n; j++)
				if (j != i && !seen[j] && sys.distsq(i,j) < sszsq)
				{	a = min(i,j);
					b = i^j^a;
					m = (bruteforce_sp_lookup[a][b] != (ui)-1 ? bruteforce_sp_lookup[a][b] : bruteforce_lookup[a][b]);
					if (m != (ui)-1)
					{	printf("run %d, pair (%d,%d) [skin]: -1 != %d\n", run, a, b, m);
						printf(" %d %d\n", bruteforce_sp_lookup[a][b], bruteforce_lookup[a][b]);
						test_fail;
					}
				}
		}
	}
	test_success;
}
