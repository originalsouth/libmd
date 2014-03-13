#ifndef rtgroups_h
#include "../../rtgroups.h"
#endif

bool test_modify_interactions()
{	rseedb = 42;
	ui runs = 100, nTypes = 100, actions = 10*nTypes, run, action, n, m, i, j, k, v;
	n = nTypes;
	ui bruteforce_type_lookup[nTypes][nTypes];
	vector<set<pair<ui,ui>>> backdoor;
	set<ui> removed;
	pair<ui,ui> id;
	md<2> sys(n);
	vector<ldf> V(1);
	for (run = 0; run < runs; run++)
	{	sys.clear();
		sys.init(n);
		memset(bruteforce_type_lookup, -1, sizeof(bruteforce_type_lookup));
		backdoor.clear();
		removed.clear();
		for (i = 0; i < n; i++)
			sys.set_type(i, i);
		m = 0;
		
		// Initialize
		for (i = 0; i < nTypes; i++)
			for (j = i; j < nTypes; j++)
				if (coinflip()) // Coinflip
				{	if (m < nTypes || randnrb() % 4 < 3)
					{	// Add new interaction
					  V[0] = m;
						if (coinflip())
							sys.add_typeinteraction(i,j,0,&V); // Directly
						else
						{	k = sys.add_interaction(0,&V); // Indirectly
							sys.add_typeinteraction(i,j,k);
						}
						bruteforce_type_lookup[i][j] = bruteforce_type_lookup[j][i] = m;
						backdoor.push_back(set<pair<ui,ui>>());
						backdoor[m].insert(make_pair(i,j));
						m++;
					}
					else
					{	// Use existing interaction
						k = randnrb() % m;
						sys.add_typeinteraction(i,j,k);
						bruteforce_type_lookup[i][j] = bruteforce_type_lookup[j][i] = v = (int)(.5+sys.network.library[k].parameters[0]);
						backdoor[v].insert(make_pair(i,j));
					}
				}
		
		// Mess around
		for (action = 0; action < actions; action++)
		{	if (randnrb() % 8 < 7)
			{	// Mess around with typeinteractions
				i = randnrb() % nTypes;
				do
					j = randnrb() % nTypes;
				while (i==j);
				if (i>j)
					swap(i,j);
				id = sys.network.hash(i,j);
				if (!sys.network.lookup.count(id))
					switch (randnrb()%2)
					{	case 0 :
							V[0] = m;
							removed.erase(sys.add_typeinteraction(i,j,0,&V));
							bruteforce_type_lookup[i][j] = m;
							backdoor.push_back(set<pair<ui,ui>>());
							backdoor[m].insert(make_pair(i,j));
							m++;
							break;
						case 1 :
							k = randnrb() % sys.network.library.size();
							if (removed.count(k))
								continue;
							sys.add_typeinteraction(i,j,k);
							bruteforce_type_lookup[i][j] = bruteforce_type_lookup[j][i] = v = (int)(.5+sys.network.library[k].parameters[0]);
							backdoor[v].insert(make_pair(i,j));
							break;
					}
				else
					switch (randnrb() % 3)
					{	case 0 :
							V[0] = m;
							sys.mod_typeinteraction(i,j,0,&V);
							backdoor[bruteforce_type_lookup[i][j]].erase(make_pair(i,j));
							bruteforce_type_lookup[i][j] = m;
							backdoor.push_back(set<pair<ui,ui>>());
							backdoor[m].insert(make_pair(i,j));
							m++;
							break;
						case 1 :
							k = randnrb() % sys.network.library.size();
							if (removed.count(k))
								continue;
							sys.mod_typeinteraction(i,j,k);
							backdoor[bruteforce_type_lookup[i][j]].erase(make_pair(i,j));
							bruteforce_type_lookup[i][j] = bruteforce_type_lookup[j][i] = v = (int)(.5+sys.network.library[k].parameters[0]);
							backdoor[v].insert(make_pair(i,j));
							break;
						case 2 :
							sys.rem_typeinteraction(i,j);
							backdoor[bruteforce_type_lookup[i][j]].erase(make_pair(i,j));
							bruteforce_type_lookup[i][j] = -1;
							break;
					}
			}
			else
			{	// Mess around with interactions
				switch (randnrb() % 3)
				{	case 0 :
						V[0] = m;
						removed.erase(sys.add_interaction(0,&V));
						backdoor.push_back(set<pair<ui,ui>>());
						m++;
						break;
					case 1 :
						V[0] = m;
						k = randnrb() % sys.network.library.size();
						v = (int)(.5+sys.network.library[k].parameters[0]);
						sys.mod_interaction(k,0,&V);
						//printf(" %d %d %d %d\n", k, v, m, (int)backdoor[v].size());
						for (auto it = backdoor[v].begin(); it != backdoor[v].end(); it++)
							bruteforce_type_lookup[it->first][it->second] = m;
						backdoor.push_back(backdoor[v]);
						backdoor[v].clear();
						m++;
						//break;
						break;
					case 2 :
						if (!sys.network.library.empty())
						{	k = randnrb() % sys.network.library.size();
							v = (int)(.5+sys.network.library[k].parameters[0]);
							sys.rem_interaction(k);
							removed.insert(k);
							for (auto it = backdoor[v].begin(); it != backdoor[v].end(); it++)
								bruteforce_type_lookup[it->first][it->second] = -1;
							backdoor[v].clear();
						}
						break;
				}
			}
		}

		// Check
		for (i = 0; i < nTypes; i++)
			for (j = i+1; j < nTypes; j++)
			{	id = sys.network.hash(i,j);
				if (sys.network.lookup.count(id))
				{	if (sys.network.library[sys.network.lookup[id]].parameters[0] != (ldf)bruteforce_type_lookup[i][j])
					{	printf("run %d, pair (%d,%d): %d != %d\n", run, i, j, (int)(.5+sys.network.library[sys.network.lookup[id]].parameters[0]), bruteforce_type_lookup[i][j]);
						test_fail
					}
				}
				else if (bruteforce_type_lookup[i][j] != (ui)-1)
				{	printf("run %d, pair (%d,%d): %d != -1\n", run, i, j, bruteforce_type_lookup[i][j]);
					test_fail
				}
			}
	}
	test_success
}
