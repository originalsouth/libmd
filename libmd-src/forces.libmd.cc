#define __libmd_src_file__
#ifndef libmd_h
#include "../libmd.h"
#endif

template<ui dim> void DAMPING(ui i,std::vector<ui> &particles,std::vector<ldf> &parameters,void *sys)
{
    //!
    //! This external damping force takes the form:
    //! \f[F^{\mu}_{\text{DAMPING}}(\dot{x}^{\mu})=-\gamma \dot{x}^{\mu}\f] <br>
    //! This function depends on one parameter:
    //! <ul>
    //! <li> the damping constant \f$\gamma\f$ </li>
    //! </ul>
    //!
    (void) particles;
    ldf gamma=parameters[0];
    for(ui d=0;d<dim;d++) SYS->particles[i].F[d]-=gamma*SYS->particles[i].dx[d];
}

template<ui dim> void DISSIPATION(ui i,std::vector<ui> &particles,std::vector<ldf> &parameters,void *sys)
{
    //!
    //! This external dissipation force takes the form:
    //! \f[F^{\mu}_{\text{DISSIPATION}}(\dot{x}^{\mu})=b \dot{x}_j^{\mu} - \dot{x}_i^{\mu}\f] <br>
    //! Here the <tt>j</tt>th particle is given in the particles std::vector<br>
    //! This function depends on one parameter:
    //! <ul>
    //! <li> the damping constant \f$b\f$ </li>
    //! </ul>
    //!
    ldf b=parameters[0];
    for(auto it: particles) for(ui d=0;d<dim;d++) SYS->particles[i].F[d]+=b*(SYS->dv(d,i,it));
}

template<ui dim> void VICSEK(ui i,std::vector<ui> &particles,std::vector<ldf> &parameters,void *sys)
{
    //!
    //! This external function implements the Vicsek model.
    //! This function depends on one parameter:
    //! <ul>
    //! <li> The cuttof radius \f$R\f$ </li>
    //! <li> The v0 the rest velocity of an active particle \f$v_{0}\f$ </li>
    //! </ul>
    //!
    static std::vector<std::vector<ldf>> dxnorm;
    if(!i)
    {
        std::vector<ldf> zero(2,0.0);
        dxnorm.assign(SYS->N,zero);
        for(ui n=0;n<SYS->N;n++)
        {
            ldf nrm=0.0;
            for(ui d=0;d<dim;d++) nrm+=pow(SYS->particles[n].dx[d],2);
            nrm=std::sqrt(nrm);
            if(nrm<std::numeric_limits<ldf>::epsilon()) for(ui d=0;d<2;d++) dxnorm[n][d]=0.0;
            else for(ui d=0;d<2;d++) dxnorm[n][d]=SYS->particles[n].dx[d]/nrm;
        }
    }
    ldf nrm=0.0;
    ldf ddx[dim];
    for(ui d=0;d<dim;d++) ddx[d]+=dxnorm[i][d];
    for(auto sij: SYS->network.skins[i]) if(SYS->distsq(i,sij.neighbor)<pow(parameters[0],2)) for(ui d=0;d<dim;d++) ddx[d]+=dxnorm[sij.neighbor][d];
    for(ui d=0;d<dim;d++) nrm+=pow(ddx[d],2);
    nrm=std::sqrt(nrm);
    for(ui d=0;d<dim;d++) SYS->particles[i].dx[d]=parameters[0]*ddx[d]/nrm;
}
