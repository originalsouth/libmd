#define __libmd_src_file__
#ifndef libmd_h
#include "../libmd.h"
#endif

template<ui dim> ldf mpmd<dim>::embedded_distsq(ui p1,ui p2)
{
    //!
    //! This function calculates the embedded distance squared between particles <tt>p1</tt> and  <tt>p2</tt>.
    //!
    using namespace std;
    return distsq(p1,p2)+pow(patch.geometryx[p2].x-patch.geometryx[p1].x,2);
}

template<ui dim> ldf mpmd<dim>::embedded_distsq(ldf x1[dim],ldf x2[dim])
{
    //!
    //! This function calculates the embedded distance squared between coordinates <tt>x1</tt> and <tt>x2</tt>.
    //!
    using namespace std;
    return distsq(x1,x2)+pow(patch.f(x2).x-patch.f(x1),2);
}

template<ui dim> ldf mpmd<dim>::embedded_distsq(ui p1,ldf x2[dim])
{
    //!
    //! This function calculates the embedded distance squared between coordinate <tt>x1</tt> and particle <tt>p2</tt>.
    //!
    using namespace std;
    return distsq(p1,x2)+pow(patch.f(x2)-patch.geometryx[p1].x,2);
}

template<ui dim> ldf mpmd<dim>::embedded_distsq(ldf x1[dim],ui p2)
{
    //!
    //! This function calculates the embedded distance squared between particle <tt>p1</tt> and coordinate <tt>x2</tt>.
    //!
    using namespace std;
    return distsq(x1,p2)+pow(patch.geometryx[p2].x-patch.f(x1),2);
}

template<ui dim> ldf mpmd<dim>::embedded_dd_p1(ui d,ui p1,ui p2)
{
    //!
    //! This function calculates the <tt>d</tt>'th component of the embedded distance vector with respect to particle p1. <br>
    //! Note this function is signed hence permutations matter.
    //!
    return dd(d,p1,p2)+((patch.geometryx[p2].x-patch.geometryx[p1].x)*patch.geometryx[p1].dx[d]);
}

template<ui dim> ldf mpmd<dim>::embedded_dd_p2(ui d,ui p1,ui p2)
{
    //!
    //! This function calculates the <tt>d</tt>'th component of the embedded distance vector with respect to particle p2. <br>
    //! Note this function is signed hence permutations matter.
    //!
    return dd(d,p2,p1)+((patch.geometryx[p1].x-patch.geometryx[p2].x)*patch.geometryx[p2].dx[d]);
}

template<ui dim> void mpmd<dim>::zuiden_C(ui i,ldf ZC[dim])
{
    //!
    //! This function calculates \f$g^{\rho \sigma} C_{\sigma} = g^{\rho \sigma} (g^p_{\sigma \mu}(x^{\mu}-x^{\mu}_p)+h^2 \frac{F^{\sigma}}{m}) \f$ for particle i of the van Zuiden integrator for particle <tt>i</tt>.
    //! It overwrites <tt>ldf ZC[dim]</tt> after it is read as input.
    //!
    using namespace std;
    ldf C[dim]={};
    memset(ZC,0,dim*sizeof(ldf));
    for(ui sigma=0;sigma<dim;sigma++) for(ui mu=0;mu<dim;mu++) C[sigma]+=patch.gp(i,sigma,mu)*(particles[i].x[mu]-particles[i].xp[mu]);
    for(ui sigma=0;sigma<dim;sigma++) C[sigma]+=pow(integrator.h,2)*particles[i].F[sigma]/particles[i].m;
    for(ui rho=0;rho<dim;rho++) for(ui sigma=0;sigma<dim;sigma++) ZC[rho]+=patch.ginv(i,rho,sigma)*C[sigma];
}

template<ui dim> void mpmd<dim>::zuiden_A(ui i,ldf eps[dim])
{
    //!
    //! This function calculates \f$g^{\rho \sigma} A_{\sigma \mu \nu} \epsilon^{\mu} \epsilon^{\nu}\f$ for particle i of the van Zuiden integrator for particle <tt>i</tt>.
    //! It overwrites <tt>ldf eps[dim]</tt> after it is read as input.
    //!
    ldf ZA[dim]={};
    for(ui rho=0;rho<dim;rho++) for(ui sigma=0;sigma<dim;sigma++) for(ui mu=0;mu<dim;mu++) for(ui nu=0;nu<dim;nu++) ZA[rho]+=patch.ginv(i,rho,sigma)*patch.A(i,sigma,mu,nu)*eps[mu]*eps[nu];
    memcpy(eps,ZA,dim*sizeof(ldf));
}

template<ui dim> void mpmd<dim>::thread_zuiden_wfi(ui i)
{
    //!
    //! This function runs the van Zuiden integrator without fixed point iterations and updates the position and velocity for particle <tt>i</tt> accordingly.
    //!
    ldf eps[dim]={};
    zuiden_C(i,eps);
    memcpy(particles[i].xp,particles[i].x,dim*sizeof(ldf));
    for(ui d=0;d<dim;d++) particles[i].x[d]+=eps[d];
    for(ui d=0;d<dim;d++) particles[i].dx[d]=eps[d]/integrator.h;
}

template<ui dim> void mpmd<dim>::thread_zuiden_protect(ui i)
{
    //!
    //! This function runs the van Zuiden integrator with fixed point iterations and updates the position and velocity for particle <tt>i</tt> accordingly.
    //! Should the integrator not converge after <tt>integrator.generations</tt> iterations this function stops with the fixed point iterations.
    //!
    using namespace std;
    ui counter=0;
    ldf ZC[dim],eps[dim],epsp[dim],val;
    zuiden_C(i,ZC);
    memcpy(eps,ZC,dim*sizeof(ldf));
    do
    {
        val=0.0;
        counter++;
        memcpy(epsp,eps,dim*sizeof(ldf));
        zuiden_A(i,eps);
        for(ui d=0;d<dim;d++) eps[d]+=ZC[d];
        for(ui d=0;d<dim;d++) val+=abs(epsp[d]-eps[d]);
        DEBUG_3("fixed point iterators cycle: " F_UI "",counter);
    }
    while(counter<integrator.generations and val>std::numeric_limits<ldf>::epsilon());
    memcpy(particles[i].xp,particles[i].x,dim*sizeof(ldf));
    for(ui d=0;d<dim;d++) particles[i].x[d]+=eps[d];
    for(ui d=0;d<dim;d++) particles[i].dx[d]=eps[d]/integrator.h;
}

template<ui dim> void mpmd<dim>::thread_zuiden(ui i)
{
    //!
    //! This function runs the van Zuiden integrator with fixed point iterations and updates the position and velocity for particle <tt>i</tt> accordingly.
    //! This function unlike mpmd<dim>::thread_zuiden_protect(ui i) waits until the integrator has converged.
    //!
    using namespace std;
    ldf ZC[dim],eps[dim],epsp[dim],val;
    zuiden_C(i,ZC);
    memcpy(eps,ZC,dim*sizeof(ldf));
    do
    {
        val=0.0;
        memcpy(epsp,eps,dim*sizeof(ldf));
        zuiden_A(i,eps);
        for(ui d=0;d<dim;d++) eps[d]+=ZC[d];
        for(ui d=0;d<dim;d++) val+=abs(epsp[d]-eps[d]);
    }
    while(val>std::numeric_limits<ldf>::epsilon());
    memcpy(particles[i].xp,particles[i].x,dim*sizeof(ldf));
    for(ui d=0;d<dim;d++) particles[i].x[d]+=eps[d];
    for(ui d=0;d<dim;d++) particles[i].dx[d]=eps[d]/integrator.h;
}

template<ui dim> void mpmd<dim>::thread_history(ui i)
{
    //!
    //! This function generates the history \f$ x^{\mu}_p \f$ for particle <tt>i</tt> from its position and velocity.
    //!
    for(ui d=0;d<dim;d++) particles[i].xp[d]=particles[i].x[d]-particles[i].dx[d]*integrator.h;
}

template<ui dim> void mpmd<dim>::history()
{
    //!
    //! This function calls mpmd<dim>::thread_history(ui i) for all particles.
    //! Additionally, it calculates and updates the important geometric derivatives.
    //!
    for(ui i=0;i<N;i++) thread_history(i);
    patch.geometryx.resize(N);
    patch.geometryxp.resize(N);
    for(ui i=0;i<N;i++) patch.calc(patch.geometryxp[i],particles[i].xp);
    for(ui i=0;i<N;i++) patch.calc(i,particles[i].x);
}

template<ui dim> void mpmd<dim>::thread_calc_geometry(ui i)
{
    //!
    //! This function update and calculates the new geometric information for a (new) particle position for particle <tt>i</tt>.
    //!
    patch.calc(i,particles[i].x);
}

template<ui dim> void mpmd<dim>::calc_geometry()
{
    //!
    //! This function calls mpmd<dim>::thread_calc_geometry(ui i) for all particles.
    //!
    for(ui i=0;i<N;i++) thread_calc_geometry(i);
}

template<ui dim> void mpmd<dim>::mp_thread_calc_pot_forces(ui i)
{
    //!
    //! This function is the Monge patch analog to md<dim>::thread_calc_pot_forces(ui i) and calculates the forces induced by the potentials acting on particle <tt>i</tt>.
    //!
    using namespace std;
    for(auto sij: network.skins[i]) if(i>sij.neighbor)
    {
        const ldf rsq=embedded_distsq(i,sij.neighbor);
        if(!network.update or rsq<pow(get_rco(sij.interaction),2))
        {
            const ldf r=sqrt(rsq);
            DEBUG_3("r = " F_LDF,r);
            const ldf dVdr=v.dr(network.library[sij.interaction].potential,r,network.library[sij.interaction].parameters);
            DEBUG_3("dV/dr = " F_LDF,dVdr);
            for(ui d=0;d<dim;d++)
            {
                particles[i].F[d]+=embedded_dd_p1(d,i,sij.neighbor)*dVdr/r;
                DEBUG_3("particles[" F_UI "].F[" F_UI "] = " F_LDF,i,d,embedded_dd_p1(d,i,sij.neighbor)*dVdr/r);
                particles[sij.neighbor].F[d]+=embedded_dd_p2(d,i,sij.neighbor)*dVdr/r;
                DEBUG_3("particles[" F_UI "].F[" F_UI "] = " F_LDF,sij.neighbor,d,embedded_dd_p2(d,i,sij.neighbor)*dVdr/r);
            }
        }
    }
}

template<ui dim> void mpmd<dim>::calc_forces()
{
    //!
    //! This function clears all the forces and then calls mpmd<dim>::mp_thread_calc_forces(ui i) for all particles.
    //!
    if(network.update and (avars.reindex or test_index()))
    {
        DEBUG_2("regenerating skinlist");
        index();
    }
    DEBUG_2("exec is here");
    avars.export_force_calc=false;
    for(ui i=0;i<N;i++) thread_clear_forces(i);
    recalc_forces();
}

template<ui dim> void mpmd<dim>::recalc_forces()
{
    //!
    //! This function calls mpmd<dim>::mp_thread_calc_forces(ui i) for all particles (without clearing the forces).
    //!
    DEBUG_3("exec is here");
    if(!network.library.empty()) for(ui i=0;i<N;i++) mp_thread_calc_pot_forces(i);
    if(!network.forcelibrary.empty()) for(ui i=0;i<N;i++) thread_calc_ext_forces(i);
}
template<ui dim> void mpmd<dim>::integrate()
{
    //!
    //! This function is the Monge patch analog to md<dim>::integrate() and calculates the particle trajectories. <br>
    //! After integrating (and updating the particle) it call md<dim>::periodicity(), and calculates the geomatric derivatives for the new particle position.
    //!
    avars.export_force_calc=true;
    switch(integrator.method)
    {
        case MP_INTEGRATOR::VVERLET:
            WARNING("flatspace integrator");
            DEBUG_2("integrating using flatspace velocity Verlet");
            for(ui i=0;i<N;i++) if(!particles[i].fix) thread_vverlet_x(i);
            recalc_forces();
            for(ui i=0;i<N;i++) if(!particles[i].fix) thread_vverlet_dx(i);
        break;
        case MP_INTEGRATOR::SEULER:
            WARNING("flatspace integrator");
            DEBUG_2("integrating using flatspace symplectic Euler");
            for(ui i=0;i<N;i++) if(!particles[i].fix) thread_seuler(i);
        break;
        case MP_INTEGRATOR::VZ_WFI:
            WARNING("incomplete integrator");
            DEBUG_2("integrating using van Zuiden without fixed point iterations");
            for(ui i=0;i<N;i++) if(!particles[i].fix) thread_zuiden_wfi(i);
        break;
        case MP_INTEGRATOR::VZ_P:
            DEBUG_2("integrating using van Zuiden with protected (with " F_UI " generations) fixed point iterations",integrator.generations);
            for(ui i=0;i<N;i++) if(!particles[i].fix) thread_zuiden_protect(i);
        break;
        default:
            DEBUG_2("integrating using van Zuiden");
            for(ui i=0;i<N;i++) if(!particles[i].fix) thread_zuiden(i);
        break;
    }
    periodicity();
    patch.geometryxp=patch.geometryx;
    calc_geometry();
}

template<ui dim> ldf mpmd<dim>::thread_T(ui i)
{
    //!
    //! This function calculates the kinetic energy of a particle <tt>i</tt> in the presence of curvature.
    //!
    ldf retval=0.0;
    for(ui mu=0;mu<dim;mu++) for(ui nu=0;nu<dim;nu++) retval+=patch.g(i,mu,nu)*particles[i].dx[mu]*particles[i].dx[nu];
    return 0.5*particles[i].m*retval;
}

template<ui dim> ldf mpmd<dim>::thread_V(ui i,bool higher_index_only)
{
    //!
    //! This function calculates the potential energy of a particle <tt>i</tt> in the presence of curvature.
    //!
    using namespace std;
    ldf retval=0.0;
    ldf rcosq=pow(network.rco,2);
    for(auto sij: network.skins[i]) if(!higher_index_only or i<sij.neighbor)
    {
        const ldf rsq=embedded_distsq(i,sij.neighbor);
        if(rsq<rcosq)
        {
            const ldf r=sqrt(rsq);
            retval+=v(network.library[sij.interaction].potential,r,network.library[sij.interaction].parameters);
            if(network.update) retval-=network.library[sij.interaction].vco;
        }
    }
    return retval;
}
