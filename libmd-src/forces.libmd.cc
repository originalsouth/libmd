#ifndef libmd_h
#include "../libmd.h"
#endif

template<ui dim> void DAMPING(ui i,vector<ui> *particles,vector<ldf> *parameters,void *sys)
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
    ldf gamma=parameters->at(0);
    for(ui d=0;d<dim;d++) SYS->particles[i].F[d]-=gamma*SYS->particles[i].dx[d];
}

template<ui dim> void DISSIPATION(ui i,vector<ui> *particles,vector<ldf> *parameters,void *sys)
{
    //!
    //! This external dissipation force takes the form:
    //! \f[F^{\mu}_{\text{DISSIPATION}}(\dot{x}^{\mu})=b \dot{x}_j^{\mu} - \dot{x}_i^{\mu}\f] <br>
    //! Here the <tt>j</tt>th particle is given in the particles vector<br>
    //! This function depends on one parameter:
    //! <ul>
    //! <li> the damping constant \f$b\f$ </li>
    //! </ul>
    //!
    ldf b=parameters->at(0);
    for(auto it: *particles) for(ui d=0;d<dim;d++) SYS->particles[i].F[d]+=b*(SYS->dv(d,i,it));
}
