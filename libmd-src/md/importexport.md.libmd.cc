#ifndef libmd_h
#include "../../libmd.h"
#endif

template<ui dim> void md<dim>::import_pos(ldf *x)
{
    //!
    //! This function takes an arbitrary number of arrays. <br>
    //! The array is assumed to be of length N, the number of particles. <br>
    //! Each <tt>i</tt>th element is then assigned to the particle the <tt>d</tt>th component of the particle's position <br>
    //! As <tt>d</tt> is a variadic variable it is stored and increased for every array.
    //!
    ui d=vvars[0];
    for(ui i=0;i<N;i++) particles[i].x[d]=x[i];
}

template<ui dim> template<typename...arg> void md<dim>::import_pos(ldf *x,arg...argv)
{
    //!
    //! This function takes an arbitrary number of arrays. <br>
    //! The array is assumed to be of length N, the number of particles. <br>
    //! Each <tt>i</tt>th element is then assigned to the particle the <tt>d</tt>th component of the particle's position <br>
    //! As <tt>d</tt> is a variadic variable it is stored and increased for every array.
    //!
    ui d=vvars[0];
    for(ui i=0;i<N;i++) particles[i].x[d]=x[i];
    import_pos(argv...);
}

template<ui dim> void md<dim>::import_pos(ui i,ldf x)
{
    //!
    //! This function takes an arbitrary number of values. <br>
    //! The values are then then assigned to the <tt>d</tt>th component of the <tt>i</tt>th particle's position <br>
    //! As <tt>d</tt> is a variadic variable it is stored and increased for every value.
    //!
    ui d=vvars[0];
    particles[i].x[d]=x;
}

template<ui dim> template<typename...arg> void md<dim>::import_pos(ui i,ldf x,arg...argv)
{
    //!
    //! This function takes an arbitrary number of values. <br>
    //! The values are then then assigned to the <tt>d</tt>th component of the <tt>i</tt>th particle's position <br>
    //! As <tt>d</tt> is a variadic variable it is stored and increased for every value.
    //!
    ui d=vvars[0];
    particles[i].x[d]=x;
    import_pos(i,argv...);
}

template<ui dim> void md<dim>::import_vel(ldf *dx)
{
    //!
    //! This function takes an arbitrary number of arrays. <br>
    //! The array is assumed to be of length N, the number of particles. <br>
    //! Each <tt>i</tt>th element is then assigned to the particle the <tt>d</tt>th component of the particle's velocity <br>
    //! As <tt>d</tt> is a variadic variable it is stored and increased for every array.
    //!
    ui d=vvars[1];
    for(ui i=0;i<N;i++) particles[i].dx[d]=dx[i];
}

template<ui dim> template<typename...arg> void md<dim>::import_vel(ldf *dx,arg...argv)
{
    //!
    //! This function takes an arbitrary number of arrays. <br>
    //! The array is assumed to be of length N, the number of particles. <br>
    //! Each <tt>i</tt>th element is then assigned to the particle the <tt>d</tt>th component of the particle's velocity <br>
    //! As <tt>d</tt> is a variadic variable it is stored and increased for every array.
    //!
    ui d=vvars[1];
    for(ui i=0;i<N;i++) particles[i].dx[d]=dx[i];
    import_vel(argv...);
}

template<ui dim> void md<dim>::import_vel(ui i,ldf dx)
{
    //!
    //! This function takes an arbitrary number of values. <br>
    //! The values are then then assigned to the <tt>d</tt>th component of the <tt>i</tt>th particle's velocity <br>
    //! As <tt>d</tt> is a variadic variable it is stored and increased for every value.
    //!
    ui d=vvars[1];
    particles[i].dx[d]=dx;
}

template<ui dim> template<typename...arg> void md<dim>::import_vel(ui i,ldf dx,arg...argv)
{
    //!
    //! This function takes an arbitrary number of values. <br>
    //! The values are then then assigned to the <tt>d</tt>th component of the <tt>i</tt>th particle's velocity<br>
    //! As <tt>d</tt> is a variadic variable it is stored and increased for every value.
    //!
    ui d=vvars[1];
    particles[i].dx[d]=dx;
    import_vel(i,argv...);
}

template<ui dim> void md<dim>::import_force(ldf *F)
{
    //!
    //! This function takes an arbitrary number of arrays. <br>
    //! The array is assumed to be of length N, the number of particles. <br>
    //! Each <tt>i</tt>th element is then assigned to the particle the <tt>d</tt>th component of the particle's force<br>
    //! As <tt>d</tt> is a variadic variable it is stored and increased for every array.
    //!
    ui d=vvars[2];
    for(ui i=0;i<N;i++) particles[i].F[d]=F[i];
}

template<ui dim> template<typename...arg> void md<dim>::import_force(ldf *F,arg...argv)
{
    //!
    //! This function takes an arbitrary number of arrays. <br>
    //! The array is assumed to be of length N, the number of particles. <br>
    //! Each <tt>i</tt>th element is then assigned to the particle the <tt>d</tt>th component of the particle's force<br>
    //! As <tt>d</tt> is a variadic variable it is stored and increased for every array.
    //!
    ui d=vvars[2];
    for(ui i=0;i<N;i++) particles[i].F[d]=F[i];
    import_force(argv...);
}

template<ui dim> void md<dim>::import_force(ui i,ldf F)
{
    //!
    //! This function takes an arbitrary number of values. <br>
    //! The values are then then assigned to the <tt>d</tt>th component of the <tt>i</tt>th particle's force <br>
    //! As <tt>d</tt> is a variadic variable it is stored and increased for every value.
    //!
    ui d=vvars[2];
    particles[i].F[d]=F;
}

template<ui dim> template<typename...arg> void md<dim>::import_force(ui i,ldf F,arg...argv)
{
    //!
    //! This function takes an arbitrary number of values. <br>
    //! The values are then then assigned to the <tt>d</tt>th component of the <tt>i</tt>th particle's force <br>
    //! As <tt>d</tt> is a variadic variable it is stored and increased for every value.
    //!
    ui d=vvars[2];
    particles[i].F[d]=F;
    import_force(i,argv...);
}

template<ui dim> void md<dim>::export_pos(ldf *x)
{
    //!
    //! This function takes an arbitrary number of arrays. <br>
    //! The array is assumed to be of length N, the number of particles. <br>
    //! Then for each <tt>i</tt>th particle position every <tt>d</tt>th component is assigned to that array <br>
    //! As <tt>d</tt> is a variadic variable it is stored and increased for every array.
    //!
    ui d=vvars[3];
    for(ui i=0;i<N;i++) x[i]=particles[i].x[d];
}

template<ui dim> template<typename...arg> void md<dim>::export_pos(ldf *x,arg...argv)
{
    //!
    //! This function takes an arbitrary number of arrays. <br>
    //! The array is assumed to be of length N, the number of particles. <br>
    //! Then for each <tt>i</tt>th particle position every <tt>d</tt>th component is assigned to that array <br>
    //! As <tt>d</tt> is a variadic variable it is stored and increased for every array.
    //!
    ui d=vvars[3];
    for(ui i=0;i<N;i++) x[i]=particles[i].x[d];
    export_pos(argv...);
}

template<ui dim> void md<dim>::export_pos(ui i,ldf &x)
{
    //!
    //! This function takes an arbitrary number of values. <br>
    //! Then for each <tt>i</tt>th particle position every <tt>d</tt>th component is assigned to that value <br>
    //! As <tt>d</tt> is a variadic variable it is stored and increased for every array.
    //!
    ui d=vvars[3];
    x=particles[i].x[d];
}

template<ui dim> template<typename...arg> void md<dim>::export_pos(ui i,ldf &x,arg...argv)
{
    //!
    //! This function takes an arbitrary number of values. <br>
    //! Then for each <tt>i</tt>th particle position every <tt>d</tt>th component is assigned to that value <br>
    //! As <tt>d</tt> is a variadic variable it is stored and increased for every array.
    //!
    ui d=vvars[3];
    x=particles[i].x[d];
    export_pos(i,argv...);
}

template<ui dim> void md<dim>::export_vel(ldf *dx)
{
    //!
    //! This function takes an arbitrary number of arrays. <br>
    //! The array is assumed to be of length N, the number of particles. <br>
    //! Then for each <tt>i</tt>th particle velocity every <tt>d</tt>th component is assigned to that array <br>
    //! As <tt>d</tt> is a variadic variable it is stored and increased for every array.
    //!
    ui d=vvars[4];
    for(ui i=0;i<N;i++) dx[i]=particles[i].dx[d];
}

template<ui dim> template<typename...arg> void md<dim>::export_vel(ldf *dx,arg...argv)
{
    //!
    //! This function takes an arbitrary number of arrays. <br>
    //! The array is assumed to be of length N, the number of particles. <br>
    //! Then for each <tt>i</tt>th particle velocity every <tt>d</tt>th component is assigned to that array <br>
    //! As <tt>d</tt> is a variadic variable it is stored and increased for every array.
    //!
    ui d=vvars[4];
    for(ui i=0;i<N;i++) dx[i]=particles[i].dx[d];
    export_vel(argv...);
}

template<ui dim> void md<dim>::export_vel(ui i,ldf &dx)
{
    //!
    //! This function takes an arbitrary number of values. <br>
    //! Then for each <tt>i</tt>th particle velocity every <tt>d</tt>th component is assigned to that value <br>
    //! As <tt>d</tt> is a variadic variable it is stored and increased for every array.
    //!
    ui d=vvars[4];
    dx=particles[i].dx[d];
}

template<ui dim> template<typename...arg> void md<dim>::export_vel(ui i,ldf &dx,arg...argv)
{
    //!
    //! This function takes an arbitrary number of values. <br>
    //! Then for each <tt>i</tt>th particle velocity every <tt>d</tt>th component is assigned to that value <br>
    //! As <tt>d</tt> is a variadic variable it is stored and increased for every array.
    //!
    ui d=vvars[4];
    dx=particles[i].dx[d];
    export_vel(i,argv...);
}

template<ui dim> void md<dim>::export_force(ldf *F)
{
    //!
    //! This function takes an arbitrary number of arrays. <br>
    //! The array is assumed to be of length N, the number of particles. <br>
    //! Then for each <tt>i</tt>th particle force every <tt>d</tt>th component is assigned to that array <br>
    //! As <tt>d</tt> is a variadic variable it is stored and increased for every array.
    //! This function checks if forces need to be recalculated
    //!
    ui d=vvars[5];
    if(avars.export_force_calc) calc_forces();
    for(ui i=0;i<N;i++) F[i]=particles[i].F[d];
}

template<ui dim> template<typename...arg> void md<dim>::export_force(ldf *F,arg...argv)
{
    //!
    //! This function takes an arbitrary number of arrays. <br>
    //! The array is assumed to be of length N, the number of particles. <br>
    //! Then for each <tt>i</tt>th particle force every <tt>d</tt>th component is assigned to that array <br>
    //! As <tt>d</tt> is a variadic variable it is stored and increased for every array.
    //! This function checks if forces need to be recalculated
    //!
    ui d=vvars[5];
    if(avars.export_force_calc) calc_forces();
    for(ui i=0;i<N;i++) F[i]=particles[i].F[d];
    export_force(argv...);
}

template<ui dim> void md<dim>::export_force(ui i,ldf &F)
{
    //!
    //! This function takes an arbitrary number of values. <br>
    //! Then for each <tt>i</tt>th particle force every <tt>d</tt>th component is assigned to that value <br>
    //! As <tt>d</tt> is a variadic variable it is stored and increased for every array.
    //!
    ui d=vvars[5];
    if(avars.export_force_calc) calc_forces();
    F=particles[i].F[d];
}

template<ui dim> template<typename...arg> void md<dim>::export_force(ui i,ldf &F,arg...argv)
{
    //!
    //! This function takes an arbitrary number of values. <br>
    //! Then for each <tt>i</tt>th particle force every <tt>d</tt>th component is assigned to that value <br>
    //! As <tt>d</tt> is a variadic variable it is stored and increased for every array.
    //!
    ui d=vvars[5];
    if(avars.export_force_calc) calc_forces();
    F=particles[i].F[d];
    export_force(i,argv...);
}

template<ui dim> ldf md<dim>::direct_readout_x(ui d,ui i)
{
    //!
    //! This function is a somewhat redundant function that returns <tt>particles[i].x[d]</tt>.
    //!
    return particles[i].x[d];
}

template<ui dim> ldf md<dim>::direct_readout_dx(ui d,ui i)
{
    //!
    //! This function is a somewhat redundant function that returns <tt>particles[i].dx[d]</tt>.
    //!
    return particles[i].dx[d];
}

template<ui dim> ldf md<dim>::direct_readout_F(ui d,ui i)
{
    //!
    //! This function returns <tt>particles[i].F[d]</tt> after it checks if a Force calculation is necessary.
    //!
    if(avars.export_force_calc) calc_forces();
    return particles[i].F[d];
}

template<ui dim> ldf md<dim>::direct_readout(ui i,uc type)
{
    //!
    //! This function is convenience function for <tt> type</tt>:
    //! <ul>
    //! <li> <tt>type='v'</tt>  md<dim>::direct_readout_dx <tt>(d=vvars,i)</tt></li>
    //! <li> <tt>type='F'</tt>  md<dim>::direct_readout_dx <tt>(d=vvars,i)</tt></li>
    //! <li> <tt>default</tt>  md<dim>::direct_readout_dx <tt>(d=vvars,i)</tt></li>
    //! </ul>
    //! vvars is variadic_var
    //!
    ui d=vvars[6];
    switch(type)
    {
        case 'v': return direct_readout_dx(d,i); break;
        case 'F': return direct_readout_F(d,i); break;
        default: return direct_readout_x(d,i); break;
    }
}

template<ui dim> ldf md<dim>::direct_readout(ui d,ui i,uc type)
{
    //!
    //! This function is convenience function for <tt> type</tt>:
    //! <ul>
    //! <li> <tt>type='v'</tt>  md<dim>::direct_readout_dx <tt>(d,i)</tt></li>
    //! <li> <tt>type='F'</tt>  md<dim>::direct_readout_dx <tt>(d,i)</tt></li>
    //! <li> <tt>default</tt>  md<dim>::direct_readout_dx <tt>(d,i)</tt></li>
    //! </ul>
    //!
    switch(type)
    {
        case 'v': return direct_readout_dx(d,i); break;
        case 'F': return direct_readout_F(d,i); break;
        default: return direct_readout_x(d,i); break;
    }
}
