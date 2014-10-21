template<ui dim> variadic_vars<dim>::variadic_vars()
{
    //!
    //! Variadic helper variables struct initializer. <br>
    //! Seven variables are initialized to zero.
    //!
    vvars.assign(7,0);
}

template<ui dim> ui variadic_vars<dim>::operator[](ui i)
{
    //!
    //! The variadic variables change their value everytime they're being read <tt>var=(var+1)%dim</tt>.
    //!
    ui retval=vvars[i];
    vvars[i]=(vvars[i]+1)%dim;
    return retval;
}

template<ui dim> void variadic_vars<dim>::reset()
{
    //!
    //! Reinitizialize the variadic variables.
    //!
    vvars.assign(7,0);
}

template<ui dim> void variadic_vars<dim>::reset(ui i)
{
    //!
    //! Reinitizialize variadic variable <tt>i</tt>.
    //!
    vvars[i]=0;
}

