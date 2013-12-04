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
