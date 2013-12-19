template<ui dim> box<dim>::box()
{
    for(ui d=0;d<dim;d++) bcond[d]=0;
    for(ui d=0;d<dim;d++) L[d]=10.0;
    for(ui d=0;d<dim;d++) { vshear[d]=0.0; xshear[d]=0.0; }
}
