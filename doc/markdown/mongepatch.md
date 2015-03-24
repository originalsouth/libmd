Monge patch
============

Curvature in libmd               {#mp-curvature}
=====================
libmd supports Molecular dynamics on manifolds (curved surfaces) in any dimension.
Although the [integrator for curved surfaces](@ref mp-integrator) supports any manifold, a Monge gauge/patch is chosen.
A Monge patch basically assumes that the curvature is introduced by using a scalar function \f$f:\mathbb{R}^{d}\rightarrow\mathbb{R}\f$, for any dimension \f$d\f$.
Such that the metric yields: 
\f{align}{
    g_{\mu \nu}= \delta_{\mu \nu} + \partial_{\mu} \partial_{\nu} f
\f}
In libmd \f$f\f$ can be controlled by the user.
Using automatic differentiation libmd calculates all the required differential geometric forms like the metric determinant, the metric inverse, and the Christoffel symbols.

STUB
STUB
STUB
