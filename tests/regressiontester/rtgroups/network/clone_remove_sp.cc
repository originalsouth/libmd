bool test_clone_remove_sp()
{	ui runs = 50, n = 50, S = 5, T = 3, nTypes = 5, nst = 5, run, mode, d, i, j, s, t;
	ldf h[2];
	ldf dr[2];
	md<2> sys(n);
	sys.simbox.L[0] = sys.simbox.L[1] = 10.0;
	sys.simbox.bcond[0] = sys.simbox.bcond[1] = BCOND::PERIODIC;
	sys.set_rco(4.0);
	sys.set_ssz(5.0);
	vector<ldf> V(2);
	for (run = 0; run < runs; run++)
	{	for (mode = 0; mode < 2; mode++)
		{	rseed = rseedb = run+42;
			sys.clear();
			sys.init(n);
			for (i = 0; i < nTypes; i++)
				for (j = i; j < nTypes; j++)
					if (randnrb() % 8 < 7)
					{	V[0] = randnr();
						V[1] = randnr();
						sys.add_typeinteraction(i,j,POT::HOOKEAN,&V);
					}
			for (t = 0; t < T; t++)
			{	sys.add_sptype();
				for (i = 0; i < nst; i++)
					for (j = i+1; j < nst; j++)
					{	V[0] = randnr();
						sys.add_sp_interaction(t,i,j,0,&V);
					}
			}
			for (s = 0; s < S; s++)
			{	sys.add_sp(randnrb() % T);
				for (j = 0; j < nst; j++)
					if (randnrb() % 8 < 7)
					{	do
							i = randnrb() % n;
						while (sys.network.spid[i] < UI_MAX);
						sys.sp_ingest(s,i,j);
					}
			}
			for (i = 0; i < n; i++)
			{	sys.set_type(i, randnrb() % nTypes);
				for (d = 0; d < 2; d++)
				{	sys.particles[i].x[d] = 10*randnr()-5.0;
					sys.particles[i].dx[d] = randnr()-0.5;
				}
			}

			// Mess around
			for (t = 0; t < 10; t++)
			{	sys.timesteps(100);
				s = randnrb() % S;
				dr[0] = randnr();
				dr[1] = randnr();
				if (mode == 1)
				{	sys.clone_sp(s, dr);
					sys.rem_sp_particles(s);
					sys.rem_sp(s);
					sys.index();
				}
				else
					sys.translate_sp(s, dr);
			}

			// Hash configuration
			h[mode] = 0;
			for (i = 0; i < sys.N; i++)
				for (d = 0; d < 2; d++)
					h[mode] += pow(sys.particles[i].x[d],3);
		}
		// Check
		if (!skins_consistent(sys.network.skins))
		{	printf("run %d: inconsistency in skins\n", run);
			test_fail;
		}
		if (!(fabs(h[0]-h[1]) <= 1e-8))
			test_fail;
	}
	test_success;
}
