integrators::integrators()
{
    method=0;
    h=1e-3;
    generation=0;
    generations=1;
}

void integrators::set(uc scheme)
{
    method=scheme;
    switch(scheme)
    {
        case 1: generations=2; break;
        case 2: generations=4; break;
        default: generations=1; break;
    }
}
