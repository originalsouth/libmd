bool test_modify_bonds()
{	rseed = rseed_stdval;
	ui runs = 100, n = 100, nTypes = 30, actions = 10*n, run, action, d, m, i, j;
	ui bruteforce_lookup[n][n];
	ui bruteforce_type_lookup[nTypes][nTypes];
	pair<ui,ui> id;
	md<2> sys(n);
	sys.simbox.L[0] = sys.simbox.L[1] = 10.0;
	sys.simbox.bcond[0] = sys.simbox.bcond[1] = BCOND::PERIODIC;
	sys.set_rco(4.0);
	sys.set_ssz(4.0);
	vector<ldf> V(1);
	for (run = 0; run < runs; run++)
	{	sys.clear();
		sys.init(n);
		memset(bruteforce_type_lookup, -1, sizeof(bruteforce_type_lookup));
		for (i = 0; i < n; i++)
		{	sys.set_type(i, irand() % nTypes);
			for (d = 0; d < 2; d++)
			{	sys.particles[i].x[d] = 10*urand()-5.0;
				sys.particles[i].dx[d] = urand()-0.5;
			}
		}
		m = 1;
		// Initialize
		for (i = 0; i < nTypes; i++)
			for (j = i; j < nTypes; j++)
				if (coinflip())
				{	V[0] = m;
					sys.add_typeinteraction(i,j,0,V);
					bruteforce_type_lookup[i][j] = bruteforce_type_lookup[j][i] = m++;
				}
		for (i = 0; i < n; i++)
			for (j = i+1; j < n; j++)
				bruteforce_lookup[i][j] = bruteforce_type_lookup[sys.particles[i].type][sys.particles[j].type];
		sys.index();
		
		// Mess around
		for (action = 0; action < actions; action++)
		{	i = irand() % n;
			do
				j = irand() % n;
			while (i==j);
			if (i>j)
				swap(i,j);
			id = sys.network.hash(sys.particles[i].type, sys.particles[j].type);
			if (!sys.network.lookup.count(id))
			{	V[0] = m;
				sys.add_bond(i,j,0,V);
				bruteforce_lookup[i][j] = m++;
			}
			else if (coinflip())
			{ V[0] = m;
				sys.mod_bond(i,j,0,V);
				bruteforce_lookup[i][j] = m++;
			}
			else
			{	sys.rem_bond(i,j);
				bruteforce_lookup[i][j] = -1;
			}
		}

		// Check
		if (!skins_consistent(sys.network.skins))
		{	printf("run %d: inconsistency in skins\n", run);
			test_fail;
		}
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
			}
	}
	test_success;
}
