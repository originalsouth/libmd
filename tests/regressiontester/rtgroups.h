#ifndef rtgroups_h
#define rtgroups_h

unsigned long long rseedb;

ui randnrb()
{   return rseedb = (16807 * rseedb) % 2147483647;
}

long long rseed;

ldf randnr()
{   return ((rseed = (16807 * rseed) % 2147483647) + .5) / 2147483647.0;
}

bool coinflip()
{   return randnrb() & 32;
}

#endif
