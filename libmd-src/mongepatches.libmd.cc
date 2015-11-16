#define __libmd_src_file__
#ifndef libmd_h
#include "../libmd.h"
#endif

ldf kdelta(ui i,ui j)
{
    //!
    //! Kronecker delta function:
    //! \f{align}{\delta_{ij}=\begin{cases} 1 &\quad i=j\\0&\quad i \neq j\end{cases}\f}
    //!
    return (i==j)?1.0:0.0;
}

template<class X,ui dim> X FLATSPACE(X x[dim],std::vector<ldf> &param)
{
    //!
    //! The trivial Monge patch:
    //! \f[f(x^{\rho})_{\text{FLATSPACE}} = 0 \f] <br>
    //! This function disregards all parameters.
    //!
    (void) x;
    (void) param;
    return 0.0;
}

template<class X,ui dim> X GAUSSIANBUMP(X x[dim],std::vector<ldf> &param)
{
    //!
    //! The Gaussian bump Monge patch:
    //! \f[f(x^{\rho})_{\text{GAUSSIANBUMP}} = A e^{-K x^{\rho}x^{\rho}}\f] <br>
    //! This function depends on two parameters:
    //! <ul>
    //! <li> the bump amplitude \f$A\f$ </li>
    //! <li> the bump width \f$K\f$ </li
    //! </ul>
    //!
    using std::pow;
    const ldf A=param[0];
    const ldf K=param[1];
    X retval=0.0;
    for(ui d=0;d<dim;d++) retval+=pow(x[d],2);
    return A*exp(-K*retval);
}

template<class X,ui dim> X EGGCARTON(X *x,std::vector<ldf> &param)
{
    //!
    //! The egg carton Monge patch:
    //! \f[f(x^{\rho})_{\text{EGGCARTON}} = A \prod^{d}_{\rho=1} \cos \left( K^{\rho} x^{\rho} \right) \f] <br>
    //! This function depends on \f$d+1\f$ parameters:
    //! <ul>
    //! <li> the amplitude \f$A\f$ </li>
    //! <li> the \f$d\f$-dimensional wave vector \f$K^{\rho}\f$ </li>
    //! </ul>
    //!
    using std::cos;
    X retval=param[0];
    for(ui d=0;d<dim;d++) retval*=cos(x[d]*param[d+1]);
    return retval;
}

template<class X,ui dim> X MOLLIFIER(X *x,std::vector<ldf> &param)
{
    //!
    //! The the mollifier Monge patch:
    //! \f{align}{f(x^{\rho})_{\text{MOLLIFIER}} = \begin{cases} Ae^{\frac{x^{\rho}x^{\rho}}{x^{\rho}x^{\rho}-K^2}} &\quad \lvert x^{\rho} \rvert < K \\0&\quad \lvert x^{\rho} \rvert \geq K \end{cases} \f} <br>
    //! This function depends on two parameters:
    //! <ul>
    //! <li> the amplitude \f$A\f$ </li>
    //! <li> the width \f$K\f$ </li>
    //! </ul>
    //!
    const ldf A=param[0];
    const ldf Ksq=std::pow(param[1],2);
    using std::pow;
    using std::exp;
    X retval=0.0;
    for(ui d=0;d<dim;d++) retval+=pow(x[d],2);
    return (retval<Ksq)?A*exp(retval/(retval-Ksq)):static_cast<X>(0.0);
}
