/* This module runs a bunch of simulations, each one twice.
 * The second time around it occasionally clones a random superparticle and removes the original one.
 * It is checked that there are no differences in the final configuration between the two runs.
 */

bool cmp (particle<2> P, particle<2> Q)
{	return P.x[0] < Q.x[0] || (P.x[0] == Q.x[0] && P.x[1] < Q.x[1]);
}

bool test_clone_remove_sp()
{	ui runs = 50, n = 50, S = 5, T = 3, nTypes = 5, nst = 5, run, mode, d, i, j, s, t;
	ldf dr[2];
	md<2> sys[2];
	for (mode = 0; mode < 2; mode++)
	{	sys[mode].simbox.L[0] = sys[mode].simbox.L[1] = 10.0;
		sys[mode].simbox.bcond[0] = sys[mode].simbox.bcond[1] = BCOND::PERIODIC;
		sys[mode].set_rco(4.0);
		sys[mode].set_ssz(5.0);
	}
	vector<ldf> V(2);
	for (run = 0; run < runs; run++)
	{	for (mode = 0; mode < 2; mode++)
		{	rseed = rseedb = run+42;
			sys[mode].clear();
			sys[mode].init(n);
			for (i = 0; i < nTypes; i++)
				for (j = i; j < nTypes; j++)
					if (randnrb() % 8 < 7)
					{	V[0] = randnr();
						V[1] = randnr();
						sys[mode].add_typeinteraction(i,j,POT::HOOKEAN,&V);
					}
			for (t = 0; t < T; t++)
			{	sys[mode].add_sptype();
				for (i = 0; i < nst; i++)
					for (j = i+1; j < nst; j++)
					{	V[0] = randnr();
						sys[mode].add_sp_interaction(t,i,j,0,&V);
					}
			}
			for (s = 0; s < S; s++)
			{	sys[mode].add_sp(randnrb() % T);
				for (j = 0; j < nst; j++)
					if (randnrb() % 8 < 7)
					{	do
							i = randnrb() % n;
						while (sys[mode].network.spid[i] < UI_MAX);
						sys[mode].sp_ingest(s,i,j);
					}
			}
			for (i = 0; i < n; i++)
			{	sys[mode].set_type(i, randnrb() % nTypes);
				for (d = 0; d < 2; d++)
				{	sys[mode].particles[i].x[d] = 10*randnr()-5.0;
					sys[mode].particles[i].dx[d] = randnr()-0.5;
				}
			}

			// Mess around
			for (t = 0; t < 10; t++)
			{	sys[mode].timesteps(100);
				s = randnrb() % S;
				dr[0] = randnr();
				dr[1] = randnr();
				if (mode == 1)
				{	sys[mode].clone_sp(s, dr);
					sys[mode].rem_sp_particles(s);
					sys[mode].rem_sp(s);
					sys[mode].index();
				}
				else
					sys[mode].translate_sp(s, dr);
			}

			// Check skins
			if (!skins_consistent(sys[mode].network.skins))
			{	printf("run %d: inconsistency in skins\n", run);
				test_fail;
			}
		}

		// Check for differences
		if (sys[0].N != sys[1].N)
			test_fail;
		for (mode = 0; mode < 2; mode++)
			sort(sys[mode].particles.begin(), sys[mode].particles.end(), cmp);
		for (i = 0; i < sys[0].N; i++)
			for (d = 0; d < 2; d++)
				if (sys[0].particles[i].type != sys[1].particles[i].type ||
				    fabs(sys[0].particles[i].x[d] - sys[1].particles[i].x[d]) > 1e-8 ||
				    fabs(sys[0].particles[i].dx[d] - sys[1].particles[i].dx[d]) > 1e-8)
					test_fail;
	}
	test_success;
}
