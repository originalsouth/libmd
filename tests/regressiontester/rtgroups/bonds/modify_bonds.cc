unsigned long long rseedb = 42;

ui randnrb()
{	return rseedb = (16807 * rseedb) % 2147483647;
}

void showall (md<2>& sys)
{	printf("\n");
	ui n = sys.N, i, j;
	pair<ui,ui> id;
	printf("         ");
	for (j = 0; j < n; j++)
		printf("%4d", j);
	printf("\n");
	for (i = 0; i < n; i++)
	{	printf("%2d (%2d)  ", i, sys.particles[i].type);
		for (j = 0; j < n; j++)
		{	id = sys.network.hash(sys.particles[i].type, sys.particles[j].type);
			if (sys.network.lookup.count(id))
				printf("%4d", (int)(.5+sys.network.library[sys.network.lookup[id]].parameters[0]));
			else
				printf("    ");
		}
		printf("\n");
	}
	printf("\n");
}

bool test_modify_bonds()
{	printf("%s: %s: ",__FILE__,__FUNCTION__);
	ui runs = 1, n = 10, t = 10, m, run, i, j, k;
	ui bruteforce_lookup[n][n];
	ui bruteforce_type_lookup[t][t];
	pair<ui,ui> id;
	md<2> sys(n);
	vector<ldf> V(1);
	for (run = 0; run < runs; run++)
	{	sys.clear();
		sys.init(n);
		memset(bruteforce_type_lookup, -1, sizeof(bruteforce_type_lookup));
		for (i = 0; i < n; i++)
			sys.set_type(i, randnrb() % t);
		m = 1;
		for (i = 0; i < t; i++)
			for (j = i; j < t; j++)
				if (randnrb() & 32) // Coinflip
				{	V[0] = m;
					sys.add_typeinteraction(i,j,0,&V);
					bruteforce_type_lookup[i][j] = bruteforce_type_lookup[j][i] = m++;
				}
		for (i = 0; i < n; i++)
			for (j = i+1; j < n; j++)
				bruteforce_lookup[i][j] = bruteforce_type_lookup[sys.particles[i].type][sys.particles[j].type];
		//showall(sys);
		for (k = 0; k < n; k++)
		{	i = randnrb() % n;
			do
				j = randnrb() % n;
			while (i==j);
			if (i>j)
				swap(i,j);
			id = sys.network.hash(sys.particles[i].type, sys.particles[j].type);
			if(!sys.network.lookup.count(id))
			{	V[0] = m;
				sys.add_bond(i,j,0,&V);
				bruteforce_lookup[i][j] = m++;
			}
			else
			{	sys.rem_bond(i,j);
				bruteforce_lookup[i][j] = -1;
			}
			//printf(" %d %d\n", i, j);
		  //showall(sys);
		  //break;
		}
		for (i = 0; i < n; i++)
			for (j = i+1; j < n; j++)
			{	id = sys.network.hash(sys.particles[i].type, sys.particles[j].type);
				if(sys.network.lookup.count(id))
				{	if (sys.network.library[sys.network.lookup[id]].parameters[0] != (ldf)bruteforce_lookup[i][j])
						return false;
				}
				else if (bruteforce_lookup[i][j] != (ui)-1)
					return false;
			}
	}
	return true;
}
