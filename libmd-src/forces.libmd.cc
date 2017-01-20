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

template<ui dim> void LANGEVIN(ui i,std::vector<ui> &particles,std::vector<ldf> &parameters,void *sys)
{
    //!
    //! This external Langevin force takes the form:
    //! \f[F^{\mu}_{\text{LANGEVIN}}(\dot{x}^{\mu})=\sqrt{\left(2 \gamma k_B T\right)} \hat{\xi}(t)\f] <br>
    //! This function depends on one parameter:
    //! <ul>
    //! <li> the temperature \f$k_B T\f$ </li>
    //! <li> the damping constant \f$\gamma\f$ </li>
    //! </ul>
    //!
    (void) particles;
    static std::random_device rd;
    static std::mt19937 mt(rd());
    static std::normal_distribution<ldf> normal(0.0,1.0);
    const ldf KbT=parameters[0];
    const ldf gamma=parameters[1];
    const ldf factor=sqrt(2.0*gamma*KbT/SYS->integrator.h);
    for(ui d=0;d<dim;d++) SYS->particles[i].F[d]+=factor*normal(mt);
}
