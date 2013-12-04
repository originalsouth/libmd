template<ui dim> particle<dim>::particle(ldf mass,ui ptype,bool fixed)
{
    m=mass;
    type=ptype;
    fix=fixed;
}
