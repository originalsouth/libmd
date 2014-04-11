bool test_remove_particles()
{	ui runs = 100, n = 20, nTypes = 5, run, mode, d, i, j, t;
	ui dim = 2;
	ldf h[2];
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
			for (i = 0; i < n; i++)
			{	sys.set_type(i, randnrb() % nTypes);
				for (d = 0; d < 2; d++)
				{	sys.particles[i].x[d] = 10*randnr()-5.0;
					sys.particles[i].dx[d] = randnr()-0.5;
				}
			}
			// Initialize
			for (i = 0; i < nTypes; i++)
				for (j = i; j < nTypes; j++)
					if (randnrb() % 8 < 7)
					{	V[0] = randnr();
						V[1] = randnr();
						sys.add_typeinteraction(i,j,POT::HOOKEAN,&V);
					}

			// Mess around
			for (t = 0; t < 10; t++)
			{	sys.timesteps(100);
				if (mode == 1)
				{	i = randnrb() % n;
					// Prepare to remove particle
					ldf m, x[dim], dx[dim];
					m = sys.particles[i].m;
					memcpy(x, sys.particles[i].x, dim*sizeof(ldf));
					memcpy(dx, sys.particles[i].dx, dim*sizeof(ldf));
					ui type = sys.particles[i].type;
					sys.rem_particle(i);
					// Re-add particle
					i = sys.add_particle(x, dx, m, type);
					sys.index();
				}
			}

			// Hash configuration
			h[mode] = 0;
			for (i = 0; i < n; i++)
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
