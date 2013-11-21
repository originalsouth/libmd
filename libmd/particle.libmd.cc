template<ui dim> particle<dim>::particle()
{
    m=1.0;
    type=0;
    fix=false;
}

template<ui dim> particle<dim>::particle(ldf mass,ui ptype,bool fixed)
{
    m=mass;
    type=ptype;
    fix=fixed;
}
