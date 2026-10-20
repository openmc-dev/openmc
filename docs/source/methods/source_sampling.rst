.. _methods_source_sampling:

===============
Source Sampling
===============

A source site specifies a particle's position :math:`\mathbf{r}`, direction
:math:`\mathbf{\Omega}`, energy :math:`E`, emission time :math:`t`, particle
type, and statistical weight :math:`w`. In a fixed source calculation, OpenMC
draws sites from the external source throughout the calculation. In an
eigenvalue calculation, the external source supplies the initial population;
subsequent generations use the fission bank as described in
:ref:`method-successive-generations`.

This chapter describes the algorithms represented by the subclasses of
:class:`openmc.SourceBase`. Input examples are given in
:ref:`usersguide_source`. Unless otherwise stated, :math:`\xi` denotes an
independent uniform random number on :math:`[0,1)`, generated as described in
:ref:`methods_random_numbers`. Probability densities below are normalized with
respect to the coordinates explicitly indicated.

.. _methods_source_selection:

----------------------------
Selecting Among Source Types
----------------------------

Suppose there are :math:`N` external sources with nonnegative strengths
:math:`s_i`, with :math:`\sum_i s_i>0`, and normalized distributions
:math:`p_i`. Their mixture is

.. math::
    :label: source-mixture

    p(\mathbf{r},\mathbf{\Omega},E,t)
    = \sum_{i=1}^{N} P_i p_i(\mathbf{r},\mathbf{\Omega},E,t),
    \qquad P_i = \frac{s_i}{\sum_j s_j}.

OpenMC first selects a source index with probability :math:`P_i` and then calls
that source's sampling method. Discrete selections of this kind use alias
tables: after initialization, each selection requires a fixed number of random
draws and table lookups, independent of the number of entries.

With :attr:`~openmc.Settings.uniform_source_sampling` enabled, the source index
is instead selected with probability :math:`1/N`. The sampled site's weight is
multiplied by

.. math::

    \frac{P_i}{1/N} = \frac{N s_i}{\sum_j s_j}.

This preserves the expected contribution of each source while allocating the
same expected number of histories to every source. These factors multiply any
weights supplied by the selected source. Source strengths determine the
relative probabilities; converting results per source particle to an absolute
physical rate is a separate normalization described in
:ref:`usersguide_tally_normalization`.

-------------------------
Independent Distributions
-------------------------

An :class:`openmc.IndependentSource` defines separate distributions for
position, direction, energy, and time. Before applying constraints, its joint
distribution factors as

.. math::

    p(\mathbf{r},\mathbf{\Omega},E,t)
    = p_r(\mathbf{r}) p_\Omega(\mathbf{\Omega}) p_E(E) p_t(t).

The particle type is fixed for each independent source. Each factor can itself
be continuous, discrete, or a mixture of distributions.

Sampling One-Dimensional Distributions
--------------------------------------

For a normalized density :math:`f(x)`, inverse-transform sampling draws
:math:`x=F^{-1}(\xi)`, where

.. math::

    F(x) = \int_{-\infty}^{x} f(u)\,du.

OpenMC provides analytic distributions as well as tabulated distributions
through :mod:`openmc.stats`. A tabulated distribution first selects an interval
using its cumulative probability, then inverts the CDF within that interval
according to the specified interpolation. For example, a linear-linear PDF on
:math:`[x_i,x_{i+1}]` has

.. math::

    f(x_i+y) = f_i + m_i y, \qquad
    F(x_i+y)-F(x_i) = f_i y + \frac{m_i}{2}y^2,

where :math:`m_i=(f_{i+1}-f_i)/(x_{i+1}-x_i)`. Inverting this quadratic gives
the sampled value. A histogram PDF has a linear CDF within each interval.
Discrete distributions select among their supplied values using their
probabilities; mixture distributions select a component and sample it.

Energy and time distributions use the same one-dimensional sampling
infrastructure. Analytic spectra, such as Maxwell and Watt spectra, use their
specialized sampling algorithms. A :class:`~openmc.stats.DecaySpectrum` selects
from the nuclide decay spectra supplied by the depletion chain, sets the source
strength from the total photon emission rate, and records the parent nuclide
with each sampled photon.

Spatial Coordinates and Volume Factors
--------------------------------------

A :class:`~openmc.stats.CartesianIndependent` distribution samples :math:`x`,
:math:`y`, and :math:`z` independently. A
:class:`~openmc.stats.CylindricalIndependent` distribution samples
:math:`(R,\phi,z)` and transforms to Cartesian coordinates. Its supplied radial
PDF is a density with respect to :math:`dR`. For a uniform density per unit
volume in a cylindrical shell, the volume element
:math:`dV=R\,dR\,d\phi\,dz` requires

.. math::

    p_R(R) = \frac{2R}{R_{\max}^2-R_{\min}^2}, \qquad
    R = \sqrt{R_{\min}^2 + \xi(R_{\max}^2-R_{\min}^2)}.

Similarly, :class:`~openmc.stats.SphericalIndependent` samples radius,
:math:`\mu=\cos\theta`, and azimuth. Uniform volume sampling in a spherical
shell requires :math:`p_r(r)\propto r^2`, giving

.. math::

    r = \left[r_{\min}^3 + \xi(r_{\max}^3-r_{\min}^3)\right]^{1/3}.

The coordinate distributions specify these radial factors explicitly. The
coordinate transformation itself introduces no additional sampling weight.
Points and boxes are special cases with fixed coordinates or uniform Cartesian
coordinates, respectively.

A :class:`~openmc.stats.PointCloud` is a discrete spatial distribution. Given
coordinates :math:`\mathbf{r}_j` with strengths :math:`a_j`, it selects
:math:`\mathbf{r}_j` with probability :math:`a_j/\sum_k a_k`. It places the
particle exactly at the selected point. If points approximate a continuous
emission density, their strengths must incorporate the volumes represented by
the points.

Angular Distributions
---------------------

An isotropic direction has density :math:`1/(4\pi)` per unit solid angle.
Because :math:`d\Omega=d\mu\,d\varphi`, it can be sampled using

.. math::

    \mu = 2\xi_1-1, \qquad \varphi=2\pi\xi_2,

and constructing the direction from the polar cosine and azimuth. A
:class:`~openmc.stats.PolarAzimuthal` distribution supplies separate PDFs for
:math:`\mu` and :math:`\varphi`, referred to a specified orientation, while
:class:`~openmc.stats.Monodirectional` always returns one direction. These
angles describe particle flight directions and are distinct from the spatial
angles used to locate a source site.

-------------------------
Sources Defined on Meshes
-------------------------

A mesh distribution first selects an element, then draws a position uniformly
within that element. If :math:`a_j` is the integrated strength assigned to
element :math:`j`, with volume :math:`V_j`, the spatial PDF is

.. math::

    p_r(\mathbf{r}) = \frac{a_j}{\sum_k a_k}\frac{1}{V_j},
    \qquad \mathbf{r}\in V_j.

An :class:`openmc.IndependentSource` with a
:class:`~openmc.stats.MeshSpatial` distribution uses the same direction, energy,
and time distributions for all elements. With ``volume_normalized=True``, the
supplied strengths represent densities per unit volume: OpenMC multiplies them
by element volumes before selecting an element. Otherwise, the supplied values
are integrated element strengths.

An :class:`openmc.MeshSource` assigns an independent source to each element.
It selects an element using that source's strength, replaces the source's
spatial distribution with uniform sampling inside the selected element, and
then samples the element's direction, energy, time, and particle type. This
allows those distributions to vary with position through the element index.
Its total strength is the sum of the element source strengths. Element
strengths are integrated strengths; a prescribed volumetric emission density
must be multiplied by element volume when constructing them.

---------------------------
Stored and Compiled Sources
---------------------------

An :class:`openmc.FileSource` loads source sites from a supported particle file.
If the file contains :math:`M` sites, each draw selects an index uniformly from
:math:`0,\ldots,M-1` and copies that site's position, direction, energy, time,
weight, and particle type. Sampling is with replacement, so repeated sites are
possible. Correlations stored in a site, such as those between position and
energy, are preserved. The stored weight is retained rather than used as the
probability of choosing that site.

An :class:`openmc.CompiledSource` calls a sampling method in a user-supplied
shared library. That method returns a complete source site and can implement
correlations or distributions beyond those exposed by the built-in classes.
Its sampling and statistical-weight corrections are the responsibility of the
custom implementation. Source selection and applicable source constraints
surround the custom sampling call as they do for other source types.

.. _methods_plasma_sources:

---------------------
Plasma Volume Sources
---------------------

The :class:`openmc.TokamakSource` and :class:`openmc.StellaratorSource` classes
sample neutron emission from flux surfaces. Both assume a nonnegative emission
density :math:`S` per unit physical volume that depends only on a radial flux
coordinate. This is an assumption of these source models; an equilibrium's
geometry alone does not determine its neutron emission profile.

For coordinates :math:`(u,\alpha,\phi)` with volume Jacobian
:math:`J=|\partial(x,y,z)/\partial(u,\alpha,\phi)|`, the spatial density in
those coordinates is

.. math::
    :label: plasma-source-joint

    p(u,\alpha,\phi) \propto S(u)J(u,\alpha,\phi).

Thus uniform angular coordinates generally do not give uniform emission per
unit volume. Define the differential volume

.. math::

    V'(u) = \int\!\int J(u,\alpha,\phi)\,d\alpha\,d\phi.

The joint PDF factors into a radial marginal and an angular conditional:

.. math::
    :label: plasma-source-factorization

    p_u(u) = \frac{S(u)V'(u)}{\int S(v)V'(v)\,dv}, \qquad
    p(\alpha,\phi\mid u) = \frac{J(u,\alpha,\phi)}{V'(u)}.

Sampling in this order accounts for the volume occupied by each flux interval.
The emission profile cancels from the conditional angular PDF.

Radially Varying Spectra
------------------------

Both plasma sources accept either one energy distribution or an energy
distribution at every radial grid point. In an interval with local coordinate
:math:`\eta\in[0,1]`, the conditional energy distribution is the mixture

.. math::

    p_E(E\mid u) = (1-\eta)p_{E,i}(E) + \eta p_{E,i+1}(E).

OpenMC selects distribution :math:`i+1` with probability :math:`\eta` and
distribution :math:`i` otherwise, then samples energy from the selected
distribution. This interpolates the distributions, rather than interpolating
two sampled energies. Both classes sample an isotropic flight direction and a
separate time distribution, which defaults to a point mass at zero. Importance
weights returned by the energy and time distributions multiply the site's
weight.

.. _methods_tokamak_source:

Tokamak Sources
---------------

The tokamak geometry uses the parameterization

.. math::

    \begin{aligned}
    R(r,\alpha) &= R_0+r\cos\psi+\Delta[1-(r/a)^2], \\
    Z(r,\alpha) &= Z_{\mathrm{shift}}+\kappa r\sin\alpha, \\
    \psi &= \alpha+\delta\sin\alpha,
    \end{aligned}

where :math:`a` is the minor radius, :math:`\kappa` is elongation,
:math:`\delta` controls triangularity, and :math:`\Delta` is the Shafranov
shift. The cylindrical azimuth :math:`\phi` is independent of these functions.
Let :math:`q=r/a`, :math:`\epsilon=a/R_0`, and :math:`D=\Delta/a`.
Differentiating the geometry gives

.. math::

    \begin{aligned}
    R_r Z_\alpha-R_\alpha Z_r
    &= \kappa r\,[G(\alpha)-2Dq\cos\alpha], \\
    G(\alpha)
    &= \cos(\delta\sin\alpha)
       +\delta\sin\alpha\cos\alpha\sin\psi.
    \end{aligned}

The coordinate volume element is :math:`R` times this poloidal Jacobian.
Up to constants independent of :math:`q` and :math:`\alpha`, the joint density
is therefore

.. math::
    :label: tokamak-source-joint

    p(q,\alpha) \propto S(q)q
    [1+\epsilon q\cos\psi+\epsilon D(1-q^2)]
    [G(\alpha)-2Dq\cos\alpha].

Sampling assumes a usable, consistently oriented map. The elongation
contributes a constant to the volume Jacobian and cancels from the normalized
PDF. The toroidal angle is uniform over the requested interval.

Radial Marginal
~~~~~~~~~~~~~~~

Integrating :eq:`tokamak-source-joint` over the poloidal angle yields

.. math::
    :label: tokamak-source-radial

    p_q(q) \propto S(q)q H(q), \qquad
    H(q) = 1+\epsilon D-\frac{3}{8}c_1\epsilon q-2\epsilon Dq^2,

where, in terms of Bessel functions of the first kind,

.. math::

    c_0 = J_0(\delta)+J_2(\delta), \qquad
    c_1 = \frac{J_1(2\delta)+J_3(2\delta)}{c_0}.

The emission profile is linearly interpolated on the user grid. OpenMC refines
each interval into at least eight subintervals with spacing at most
:math:`10^{-3}`, evaluates :eq:`tokamak-source-radial` there, and constructs a
linear-linear tabulated distribution. Radial sampling then uses CDF inversion
of that tabulation. For circular, concentric surfaces,
:math:`\delta=\Delta=0`, the radial marginal reduces to
:math:`p_q(q)\propto S(q)q`.

Conditional Poloidal Sampling
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The two bracketed factors in :eq:`tokamak-source-joint` are quadratic and
linear in :math:`q`. Write their Bernstein forms as

.. math::

    \begin{aligned}
    \frac{R}{R_0} &= b_0(1-q)^2+2b_1q(1-q)+b_2q^2, \\
    G-2Dq\cos\alpha &= b_3(1-q)+b_4q,
    \end{aligned}

with

.. math::

    \begin{aligned}
    b_0 &= 1+\epsilon D, & b_1 &= b_0+\tfrac12\epsilon\cos\psi, \\
    b_2 &= 1+\epsilon\cos\psi, & b_3 &= G(\alpha), \\
    b_4 &= G(\alpha)-2D\cos\alpha.
    \end{aligned}

Their product is a sum of six separable terms
:math:`\sum_{k=0}^5 w_k(q)g_k(\alpha)`, where

.. math::

    \begin{aligned}
    (g_0,g_1,g_2,g_3,g_4,g_5)
      &= (b_0b_3,b_1b_3,b_2b_3,b_0b_4,b_1b_4,b_2b_4), \\
    (w_0,w_1,w_2,w_3,w_4,w_5)
      &= ((1-q)^3,\,2q(1-q)^2,\,q^2(1-q), \\
      &\qquad q(1-q)^2,\,2q^2(1-q),\,q^3).
    \end{aligned}

Define :math:`I_k=\int_0^\pi g_k(\alpha)\,d\alpha` and scaled integrals
:math:`\widehat I_k=I_k/(\pi c_0)`. Analytic integration gives

.. math::

    \begin{aligned}
    \widehat I_0 &= 1+\epsilon D, \\
    \widehat I_1 &= 1+\epsilon D-\tfrac{3}{16}c_1\epsilon, \\
    \widehat I_2 &= 1-\tfrac38 c_1\epsilon, \\
    \widehat I_3 &= 1+\epsilon D, \\
    \widehat I_4 &= 1+\tfrac12\epsilon D-\tfrac{3}{16}c_1\epsilon, \\
    \widehat I_5 &= 1-\epsilon D-\tfrac38 c_1\epsilon.
    \end{aligned}

Consequently, the conditional PDF on :math:`[0,\pi]` is a mixture of the six
normalized densities :math:`g_k/I_k`, with component probabilities
:math:`w_k(q)I_k/\sum_j w_j(q)I_j`. The implementation tabulates each component
using ``n_alpha`` points and samples its tabulated CDF after component
selection. Up-down symmetry supplies the other half-plane by
replacing :math:`\alpha` with :math:`2\pi-\alpha` with probability one half.
The position is then mapped to :math:`(R\cos\phi,R\sin\phi,Z)`. This
construction samples the spatial coordinates without a rejection loop, with
discretization error from the radial and poloidal tabulations.

.. note::

    The continuous integrals above give the required mixture weights. The
    current implementation stores different values for :math:`I_2`,
    :math:`I_4`, and :math:`I_5`; that mismatch introduces an additional error
    in poloidal sampling beyond the tabulation error.

.. _methods_stellarator_source:

Stellarator Sources
-------------------

Fourier Representation and Volume Jacobian
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The stellarator source uses flux coordinates :math:`(\rho,\theta,\zeta)`,
where :math:`\rho=\sqrt{s}` for normalized toroidal flux :math:`s`,
:math:`\theta` is a poloidal angle, and :math:`\zeta` equals the cylindrical
azimuth. With :math:`N_{\mathrm{fp}}` field periods, define
:math:`\chi_k=m_k\theta-n_kN_{\mathrm{fp}}\zeta`. The surfaces are

.. math::
    :label: stellarator-source-fourier

    \begin{aligned}
    R(\rho,\theta,\zeta)
      &= \sum_k [R_k^c(\rho)\cos\chi_k+R_k^s(\rho)\sin\chi_k], \\
    Z(\rho,\theta,\zeta)
      &= \sum_k [Z_k^s(\rho)\sin\chi_k+Z_k^c(\rho)\cos\chi_k].
    \end{aligned}

The sine coefficients of :math:`R` and cosine coefficients of :math:`Z` vanish
for stellarator-symmetric geometry. The same sampling procedure handles both
symmetric and asymmetric tables.

The Cartesian position is
:math:`\mathbf{r}=(R\cos\zeta,R\sin\zeta,Z)`. Applying the chain rule through
cylindrical coordinates gives

.. math::
    :label: stellarator-source-jacobian

    \left|\frac{\partial(x,y,z)}{\partial(\rho,\theta,\zeta)}\right|
    = R|\tau|, \qquad
    \tau=R_\rho Z_\theta-R_\theta Z_\rho.

Here subscripts denote partial derivatives. The derivatives :math:`R_\zeta`
and :math:`Z_\zeta` drop out of the determinant because the cylindrical azimuth
is exactly :math:`\zeta`. Thus the joint coordinate PDF is
:math:`p(\rho,\theta,\zeta)\propto S(\rho)R|\tau|`.

The map is assumed to have positive :math:`R`, nested surfaces, and a
consistently oriented Jacobian away from the magnetic axis. Either orientation
is allowed. Initialization determines a global sign :math:`\sigma` from a
signed volume estimate and evaluates :math:`f=\sigma R\tau`, so that
:math:`f=R|\tau|` for a valid map. Grid checks reject nonpositive radii or
significant changes in Jacobian sign; these checks are not a global geometric
certification.

Interpolation Within a Radial Interval
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

On interval :math:`i`, define its width :math:`h_i=\rho_{i+1}-\rho_i` and local
coordinate :math:`\eta=(\rho-\rho_i)/h_i`. Each Fourier coefficient :math:`C`
and the emission profile are interpolated linearly:

.. math::

    C(\eta)=(1-\eta)C_i+\eta C_{i+1}, \qquad
    S(\eta)=S_i+\eta(S_{i+1}-S_i).

Because Fourier evaluation is linear in the coefficients, at fixed angles
:math:`R`, :math:`Z`, :math:`R_\theta`, and :math:`Z_\theta` are linear in
:math:`\eta`, while

.. math::

    R_\rho=\frac{R_{i+1}-R_i}{h_i}, \qquad
    Z_\rho=\frac{Z_{i+1}-Z_i}{h_i}

are constant within the interval. It follows that :math:`\tau` is linear in
:math:`\eta` and :math:`f=\sigma R\tau` is quadratic. The differential volume
is the angular integral of this density. Denote its value at
:math:`\rho=\rho_i+h_i\eta` by :math:`V'_i(\eta)`:

.. math::
    :label: stellarator-source-volume-polynomial

    \begin{aligned}
    V'_i(\eta)
      &\equiv \left.\frac{dV}{d\rho}\right|_{\rho=\rho_i+h_i\eta} \\
      &= \int_0^{2\pi}\!\int_0^{2\pi}
         f_i(\eta,\theta,\zeta)\,d\theta\,d\zeta \\
      &= A_0+A_1\eta+A_2\eta^2.
    \end{aligned}

Angular integration preserves the quadratic dependence on :math:`\eta`, so
the integral and the polynomial are two expressions for the same function.
The coefficients :math:`A_0`, :math:`A_1`, and :math:`A_2` are constant within
each interval. The prime denotes differentiation with respect to :math:`\rho`;
:math:`\eta` specifies where that derivative is evaluated. Since
:math:`d\rho=h_i\,d\eta`, the corresponding volume element is
:math:`dV=h_i V'_i(\eta)\,d\eta`.

The derivative generally has different one-sided values at an interior radial
surface because the coefficient slopes change there. Each interval uses its
own values of :math:`R_\rho` and :math:`Z_\rho`; the endpoint volumes from
adjacent intervals are kept separate.

Computing Differential Volumes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To determine the three polynomial coefficients in
:eq:`stellarator-source-volume-polynomial`, initialization computes the angular
integral at just three radial positions, :math:`\eta=0,1/2,1`. At each of
these positions, it evaluates :math:`f_i` on a uniform angular grid with sizes

.. math::

    N_\theta=\max(64,8m_{\max}), \qquad
    N_\zeta=\max(64,8|n|_{\max}).

The poloidal grid covers :math:`[0,2\pi)` and the toroidal grid covers one field
period, :math:`[0,2\pi/N_{\mathrm{fp}})`. The angular integral defining
:math:`V'_i` is evaluated by the quadrature sum

.. math::

    V'_i(\eta)
    = \frac{4\pi^2}{N_\theta N_\zeta}
      \sum_{j=0}^{N_\theta-1}\sum_{k=0}^{N_\zeta-1}
      f_i(\eta,\theta_j,\zeta_k).

The weight includes all field periods. Exactness of this quadrature follows
from discrete Fourier orthogonality, or the roots-of-unity identity:

.. math::

    \frac{1}{N}\sum_{j=0}^{N-1} e^{2\pi\mathrm{i}\ell j/N}
    = \begin{cases}
        1, & \ell \text{ is a multiple of } N, \\
        0, & \text{otherwise}.
      \end{cases}

This identity follows from summing a finite geometric series. For a
trigonometric polynomial with modes :math:`-K,\ldots,K`, choosing
:math:`N>K` makes every nonconstant mode sum to zero. The periodic trapezoidal
rule then returns exactly the integral of the constant mode. This result is
stated in Corollary 3.3 of [TrefethenWeideman]_.

Products in :math:`R\tau` have poloidal and toroidal degrees at most three
times the respective geometry mode limits. Applying the identity independently
in both angles gives sufficient conditions for exact integration:

.. math::

    N_\theta>3m_{\max}, \qquad N_\zeta>3|n|_{\max}.

The grid sizes above satisfy these conditions. The quadrature is thus exact
in exact arithmetic, provided the Jacobian has a consistent sign. Recovering
only the integral requires fewer points than reconstructing the full Fourier
series under the Nyquist sampling condition.

Let :math:`v_0=V'_i(0)`, :math:`v_m=V'_i(1/2)`, and :math:`v_1=V'_i(1)` be
the three values obtained by quadrature. These determine the coefficients:

.. math::

    A_0=v_0, \qquad
    A_1=4v_m-3v_0-v_1, \qquad
    A_2=2(v_0+v_1-2v_m).

The resulting polynomial supplies :math:`V'_i(\eta)` at any radius within the
interval without further angular integration.

Only two adjacent surfaces' angular grids are held in memory at once. Their
positions and poloidal derivatives give the interval's radial derivatives by
finite differences. Initialization retains the per-interval integrals and
angular bounds needed for sampling.

Selecting and Sampling a Radial Interval
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Multiplying the linear emission profile by
:eq:`stellarator-source-volume-polynomial` yields an unnormalized cubic radial
PDF in the local coordinate:

.. math::
    :label: stellarator-source-radial-polynomial

    \begin{aligned}
    P_i(\eta) &= S(\eta)V'_i(\eta)
                 = c_0+c_1\eta+c_2\eta^2+c_3\eta^3, \\
    c_0 &= S_i A_0, \\
    c_1 &= S_i A_1+\Delta S_i A_0, \\
    c_2 &= S_i A_2+\Delta S_i A_1, \\
    c_3 &= \Delta S_i A_2,
    \end{aligned}

where :math:`\Delta S_i=S_{i+1}-S_i`. Its integrated mass includes the change
of variable :math:`d\rho=h_i\,d\eta`:

.. math::
    :label: stellarator-source-interval-mass

    \begin{aligned}
    M_i &= h_i\int_0^1 P_i(\eta)\,d\eta \\
        &= h_i\left(\frac{S_i v_0}{6}
           +\frac{(S_i+S_{i+1})v_m}{3}
           +\frac{S_{i+1}v_1}{6}\right).
    \end{aligned}

The last expression is Simpson's rule, which is exact for a cubic. OpenMC
selects interval :math:`i` with probability :math:`M_i/\sum_j M_j` using an
alias table. A common scaling of all emission values cancels from these
probabilities; the implementation divides them by their maximum before forming
the masses and PDFs.

Within the selected interval, OpenMC uses uniform proposals
:math:`\eta=\xi_1` with rejection sampling. An upper bound follows by writing
the cubic in the Bernstein basis:

.. math::

    P_i(\eta) = \sum_{j=0}^{3}\beta_j
       {3\choose j}\eta^j(1-\eta)^{3-j},

where

.. math::

    \begin{aligned}
    \beta_0 &= c_0, \\
    \beta_1 &= c_0+c_1/3, \\
    \beta_2 &= c_0+2c_1/3+c_2/3, \\
    \beta_3 &= c_0+c_1+c_2+c_3.
    \end{aligned}

The Bernstein basis functions are nonnegative and sum to one, so
:math:`B_i^{r}=\max_j\beta_j` bounds :math:`P_i` throughout the interval.
The implementation adds a small relative margin for floating-point roundoff
and stores the coefficients divided by that bound. A proposal is accepted when

.. math::

    \xi_2 < \frac{P_i(\eta)}{B_i^{r}}.

On rejection, only :math:`\eta` is redrawn; the selected interval is retained.
Selecting a new interval after a rejection would change its probability in
proportion to its acceptance rate. After acceptance,
:math:`\rho=\rho_i+h_i\eta` has the required marginal
:math:`p_\rho\propto S V'`.

Conditional Angular Rejection
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

With :math:`\rho` fixed, the target PDF is

.. math::

    p(\theta,\zeta\mid\rho)
      = \frac{f_i(\eta,\theta,\zeta)}{V'_i(\eta)}.

At each angular grid point, the three radial evaluations determine the
quadratic :math:`f_i(\eta)` throughout the interval. Its maximum is obtained
from the endpoints and, for a concave quadratic, an interior vertex. The current
angular envelope is

.. math::

    B_i^{\mathrm{ang}} = 1.05\max_{j,k}\max_{0\le\eta\le1}
       f_i(\eta,\theta_j,\zeta_k).

The radial maximization here is exact in exact arithmetic. The angular maximum
is estimated on a grid, with a 5% heuristic allowance for peaks between nodes.
Exact Fourier quadrature of an integral does not establish an upper bound on
the function. This envelope therefore does not provide a rigorous bound for
every supported coefficient table; correctness of angular rejection requires
it to bound the continuous angular density.

The sampler proposes :math:`\theta=2\pi\xi_1` and
:math:`\zeta=2\pi\xi_2` over the full torus, and accepts when

.. math::

    \xi_3 B_i^{\mathrm{ang}} < f_i(\eta,\theta,\zeta).

Rejected angles are redrawn while retaining both the radial interval and
:math:`\rho`. If the envelope is valid, this produces the angular conditional
PDF above. An underestimated envelope clips the peaks because acceptance
cannot exceed one; the current sampling loop has no explicit check for that
condition. The angular rejection limit detects excessive rejection, which is a
different failure condition.

Sampling Sequence and Approximation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The steps in the stellarator source's sampling method are:

1. Select a radial interval using its integrated mass :math:`M_i`.
2. Sample the local radial coordinate by cubic-PDF rejection within that
   interval, and calculate :math:`\rho`.
3. Sample the poloidal and toroidal angles by conditional rejection at that
   fixed radius.
4. Map the accepted position to :math:`(R\cos\zeta,R\sin\zeta,Z)`.
5. Sample an isotropic flight direction.
6. Sample energy from the common spectrum or the mixture of spectra bracketing
   :math:`\rho`, then sample emission time and apply their weight factors.

Radial integration and rejection sampling are consistent with the piecewise
linear coefficient and emission tables, up to floating-point effects. Those
tables approximate the underlying equilibrium and emission profile; increasing
the number of radial surfaces can reduce interpolation error. The finite
Fourier representation also limits the geometric detail retained. These input
approximations are distinct from the angular-envelope condition above.

The :meth:`~openmc.StellaratorSource.from_vmec` and
:meth:`~openmc.StellaratorSource.from_desc` constructors prepare the coefficient
tables before C++ initialization. VMEC's surfaces, tabulated in :math:`s`, are
relabeled by :math:`\rho=\sqrt{s}`. For DESC, the Fourier-Zernike radial basis
is evaluated at the requested radial surfaces and the angular product basis is
converted to :eq:`stellarator-source-fourier`. Both constructors convert
geometry coefficients from meters to centimeters. The user supplies the
emission profile and energy distributions separately.

.. _methods_source_constraints:

-------------------------------
Constraints and Biased Sampling
-------------------------------

Source constraints restrict sites to specified cells, materials, or universes,
to fissionable material, or to specified energy and time bounds. Positions
must also lie within the model geometry. For an acceptance region
:math:`\mathcal{A}`, resampling rejected candidates produces the conditional
distribution

.. math::

    p_{\mathcal{A}}(X) =
    \frac{p(X)\,\mathbf{1}_{\mathcal{A}}(X)}
         {\int_{\mathcal{A}}p(Y)\,dY},

where :math:`X` denotes the sampled phase-space variables. The base source
wrapper repeats sampling until a candidate passes, or, with the ``kill``
strategy, returns a failing candidate with weight zero. Independent sources
handle constraints internally: they resample positions for spatial constraints
and energies for energy constraints and the available energy range. Their time
samples are drawn directly. Constraints on an element source in a
:class:`openmc.MeshSource` are applied after selecting that element; constraints
on the mesh source itself are checked by its outer wrapper. Similarly, selecting
among multiple top-level sources precedes each source's constraint handling.
Resampling therefore conditions each selected source separately and retains its
selection probability.

OpenMC limits the number of consecutive rejections and monitors rejected sites
to detect source definitions that cannot efficiently produce acceptable
particles. The geometry and source domains should be chosen to keep acceptance
reasonable.

Source biasing draws from a proposal distribution :math:`g` in place of the
desired density :math:`p`, with the correction
:math:`w\leftarrow w\,p(X)/g(X)`. The proposal must cover every region where
the desired distribution contributes. When multiple independent factors are
biased, their corrections multiply. Mesh and point-cloud distributions can
similarly change their discrete selection probabilities and apply the ratio
of the original to the biased normalized probability. The derivation and
available biasing mechanisms are described in :ref:`methods_source_biasing`.

.. [TrefethenWeideman] Lloyd N. Trefethen and J. A. C. Weideman,
   `The Exponentially Convergent Trapezoidal Rule
   <https://people.maths.ox.ac.uk/trefethen/sirev56-3_385.pdf>`_,
   *SIAM Review*, **56** (3), 385--458 (2014), Corollary 3.3, p. 398.
