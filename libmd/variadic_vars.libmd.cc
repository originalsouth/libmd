template<ui dim> variadic_vars<dim>::variadic_vars()
{
    vvars.assign(6,0);
}

template<ui dim> ui variadic_vars<dim>::operator[](ui i)
{
    ui retval=vvars[i];
    vvars[i]=(vvars[i]+1)%dim;
    return retval;
}
