threads::threads()
{
    nothreads=1;
    block.resize(nothreads);
}

threads::threads(ui nrthreads)
{
    nothreads=nrthreads;
    block.resize(nothreads);
}

void threads::set(ui nrthreads)
{
    nothreads=nrthreads;
    block.clear();
    block.resize(nothreads);
}

void threads::setmax()
{
    nothreads=thread::hardware_concurrency();
    block.clear();
    block.resize(nothreads);
}
