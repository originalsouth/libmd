/* This module runs a bunch of simulations, each one twice.
 * The second time around it occasionally removes a random particle and inserts an identical one.
 * It is checked that there are no differences in the final configuration between the two runs.
 */

bool test_remove_particles()
{	ui runs = 100, n = 20, nTypes = 5, run, mode, d, i, j, t;
	ui dim = 2;
	ldf dr[dim];
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
			for (i = 0; i < n; i++)
			{	sys[mode].set_type(i, randnrb() % nTypes);
				for (d = 0; d < 2; d++)
				{	sys[mode].particles[i].x[d] = 10*randnr()-5.0;
					sys[mode].particles[i].dx[d] = randnr()-0.5;
				}
			}
			// Initialize
			for (i = 0; i < nTypes; i++)
				for (j = i; j < nTypes; j++)
					if (randnrb() % 8 < 7)
					{	V[0] = randnr();
						V[1] = randnr();
						sys[mode].add_typeinteraction(i,j,POT::HOOKEAN,&V);
					}

			// Mess around
			for (t = 0; t < 10; t++)
			{	sys[mode].timesteps(100);
				if (run % 2 == 0 && mode == 1)
				{	i = randnrb() % n;
					// Prepare to remove particle
					ldf m, x[dim], dx[dim];
					m = sys[mode].particles[i].m;
					memcpy(x, sys[mode].particles[i].x, dim*sizeof(ldf));
					memcpy(dx, sys[mode].particles[i].dx, dim*sizeof(ldf));
					ui type = sys[mode].particles[i].type;
					sys[mode].rem_particle(i);
					// Re-add particle
					i = sys[mode].add_particle(x, dx, m, type);
					sys[mode].index();
				}
				else if (run % 2 == 1)
				{	i = randnrb() % n;
					for (d = 0; d < dim; d++)
						dr[d] = randnr();
					if (mode == 1) // Clone particle and remove original
					{	sys[mode].clone_particle(i, dr);
						sys[mode].rem_particle(i);
						sys[mode].index();
					}
					else // Translate
						sys[mode].translate_particle(i, dr);
				}
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
			sort(sys[mode].particles.begin(), sys[mode].particles.end(), compare_particles());
		for (i = 0; i < sys[0].N; i++)
			for (d = 0; d < 2; d++)
				if (sys[0].particles[i].type != sys[1].particles[i].type ||
				    fabs(sys[0].particles[i].x[d] - sys[1].particles[i].x[d]) > 1e-8 ||
				    fabs(sys[0].particles[i].dx[d] - sys[1].particles[i].dx[d]) > 1e-8)
					test_fail;
	}
	test_success;
}
