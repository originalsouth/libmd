bool test_modify_sp_interactions()
{	rseed = rseed_stdval;
	ui runs = 100, n = 100, S = 10, T = 3, nst = 10, nTypes = 30, actions = 10*n;
	ui run, action, d, s, t, m, i, j, k, v;
	ui bruteforce_sp_lookup[n][n];
	ui bruteforce_sptype_lookup[T][nst][nst];
	bool sp_pos_used[S][nst];
	set<ui> removed;
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
		memset(bruteforce_sptype_lookup, -1, sizeof(bruteforce_sptype_lookup));
		memset(bruteforce_sp_lookup, -1, sizeof(bruteforce_sp_lookup));
		memset(sp_pos_used, false, sizeof(sp_pos_used));
		removed.clear();
		
		// Initialize
		m = 1;
		for (t = 0; t < T; t++)
		{ sys.add_sptype();
			for (i = 0; i < nst; i++)
				for (j = i+1; j < nst; j++)
					if (coinflip())
					{	V[0] = m;
						sys.add_sp_interaction(t,i,j,0,V);
						bruteforce_sptype_lookup[t][i][j] = bruteforce_sptype_lookup[t][j][i] = m++;
					}
		}
		for (s = 0; s < S; s++)
			sys.add_sp(irand() % T);
		for (i = 0; i < n; i++)
		{	sys.set_type(i, irand() % nTypes);
			if (coinflip())
			{ do
				{	s = irand() % S;
					j = irand() % nst;
				}
				while (sp_pos_used[s][j]);
				sp_pos_used[s][j] = true;
				sys.sp_ingest(s,i,j);
			}
			for (d = 0; d < 2; d++)
			{	sys.particles[i].x[d] = 10*urand()-5.0;
				sys.particles[i].dx[d] = urand()-0.5;
			}
		}
		sys.index();

		
		// Mess around
		for (action = 0; action < actions; action++)
		{	// Mess around with typeinteractions
			t = irand() % T;
			i = irand() % nst;
			do
				j = irand() % nst;
			while (i==j);
			if (i>j)
				swap(i,j);
			id = sys.network.hash(i,j);
			if (!sys.network.sptypes[t].splookup.count(id))
				switch (irand() % 2)
				{	case 0 : // Add new typeinteraction
						V[0] = m;
						removed.erase(sys.add_sp_interaction(t,i,j,0,V));
						bruteforce_sptype_lookup[t][i][j] = m;
						m++;
						break;
					case 1 : // Use existing typeinteraction
						k = irand() % sys.network.library.size();
						if (removed.count(k))
							continue;
						sys.add_sp_interaction(t,i,j,k);
						bruteforce_sptype_lookup[t][i][j] = bruteforce_sptype_lookup[t][j][i] = v = (int)(.5+sys.network.library[k].parameters[0]);
						break;
				}
			else
				switch (irand() % 3)
				{	case 0 : // Modify typeinteraction to new one
						V[0] = m;
						sys.mod_sp_interaction(t,i,j,0,V);
						bruteforce_sptype_lookup[t][i][j] = m;
						m++;
						break;
					case 1 : // Modify typeinteraction to existing one
						k = irand() % sys.network.library.size();
						if (removed.count(k))
							continue;
						sys.mod_sp_interaction(t,i,j,k);
						bruteforce_sptype_lookup[t][i][j] = bruteforce_sptype_lookup[t][j][i] = v = (int)(.5+sys.network.library[k].parameters[0]);
						break;
					case 2 : // Remove typeinteraction
						sys.rem_sp_interaction(t,i,j);
						bruteforce_sptype_lookup[t][i][j] = -1;
						break;
				}
		}

		// Check
		for (t = 0; t < T; t++)
			for (i = 0; i < nst; i++)
				for (j = i+1; j < nst; j++)
				{	id = sys.network.hash(i,j);
					if (sys.network.sptypes[t].splookup.count(id))
					{	if (sys.network.library[sys.network.sptypes[t].splookup[id]].parameters[0] != (ldf)bruteforce_sptype_lookup[t][i][j])
						{	printf("run %d, sptype %d, pair (%d,%d): %d != %d\n", run, t, i, j, (int)(.5+sys.network.library[sys.network.sptypes[t].splookup[id]].parameters[0]), bruteforce_sptype_lookup[t][i][j]);
							test_fail;
						}
					}
					else if (bruteforce_sptype_lookup[t][i][j] != (ui)-1)
					{	printf("run %d, sptype %d, pair (%d,%d): -1 != %d\n", run, t, i, j, bruteforce_sptype_lookup[t][i][j]);
						test_fail;
					}
				}
	}
	test_success;
}
