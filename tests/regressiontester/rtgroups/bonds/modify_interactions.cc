#ifndef rtgroups_h
#include "../../rtgroups.h"
#endif

/*void showall (md<2>& sys)
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
}*/

bool test_modify_interactions()
{	rseedb = 42;
	ui runs = 100, n, t = 100, m, run, i, j, k;
	n = t;
	//ui bruteforce_lookup[n][n];
	ui bruteforce_type_lookup[t][t];
	pair<ui,ui> id;
	md<2> sys(n);
	vector<ldf> V(1);
	for (run = 0; run < runs; run++)
	{	sys.clear();
		sys.init(n);
		memset(bruteforce_type_lookup, -1, sizeof(bruteforce_type_lookup));
		for (i = 0; i < n; i++)
			sys.set_type(i, i);
		m = 0;
		for (i = 0; i < t; i++)
			for (j = i; j < t; j++)
				if (randnrb() & 32) // Coinflip
				{	if (m < t || randnrb() & 40)
					{	// Add new interaction
					  V[0] = m;
						if (randnrb() & 32)
							sys.add_typeinteraction(i,j,0,&V); // Directly
						else
						{	k = sys.add_interaction(0,&V); // Indirectly
							sys.add_typeinteraction(i,j,k);
						}
						bruteforce_type_lookup[i][j] = bruteforce_type_lookup[j][i] = m++;
					}
					else
					{	// Use existing interaction
						k = randnrb() % m;
						sys.add_typeinteraction(i,j,k);
						bruteforce_type_lookup[i][j] = bruteforce_type_lookup[j][i] = (int)(.5+sys.network.library[k].parameters[0]);
					}
				}
		for (k = 0; k < n; k++)
		{	i = randnrb() % t;
			do
				j = randnrb() % t;
			while (i==j);
			if (i>j)
				swap(i,j);
			id = sys.network.hash(i,j);
			if(!sys.network.lookup.count(id))
			{	V[0] = m;
				sys.add_typeinteraction(i,j,0,&V);
				bruteforce_type_lookup[i][j] = m++;
			}
			else if (randnrb() & 32)
			{	V[0] = m;
				sys.mod_typeinteraction(i,j,0,&V);
				bruteforce_type_lookup[i][j] = m++;
			}
			else
			{	sys.rem_typeinteraction(i,j);
				bruteforce_type_lookup[i][j] = -1;
			}
		}
		for (i = 0; i < t; i++)
			for (j = i+1; j < t; j++)
			{	id = sys.network.hash(i,j);
				if(sys.network.lookup.count(id))
				{	if (sys.network.library[sys.network.lookup[id]].parameters[0] != (ldf)bruteforce_type_lookup[i][j])
						test_fail
				}
				else if (bruteforce_type_lookup[i][j] != (ui)-1)
					test_fail
			}
	}
	test_success
}
