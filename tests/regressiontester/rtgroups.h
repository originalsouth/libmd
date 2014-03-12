#ifndef rtgroups_h
#define rtgroups_h

unsigned long long rseedb = 42;

ui randnrb()
{   return rseedb = (16807 * rseedb) % 2147483647;
}

long long rseed = 42;

ldf randnr()
{   return ((rseed = (16807 * rseed) % 2147483647) + .5) / 2147483647.0;
}

#endif
