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
:ref:`methods_random_numbers`. Normalized probability densities are denoted by
:math:`p` or :math:`f`, cumulative distribution functions by :math:`F`, and
discrete selection probabilities by :math:`P`. A tilde on :math:`p` indicates
an unnormalized density. Geometric volume factors are denoted by :math:`J`.
Densities are defined with respect to the coordinates explicitly indicated.

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

For a scalar random variable :math:`X` with cumulative distribution function
:math:`F(x)=\Pr(X\le x)`, inverse-transform sampling draws
:math:`x=F^{-1}(\xi)`. If the variable has a normalized density :math:`f`, then

.. math::

    F(x) = \int_{-\infty}^{x} f(u)\,du.

OpenMC provides analytic distributions as well as tabulated distributions
through :mod:`openmc.stats`. A tabulated distribution first selects an interval
using its cumulative probability, then inverts the CDF within that interval
according to the specified interpolation. For example, a linear-linear PDF on
:math:`[x_i,x_{i+1}]`, with :math:`f_i=f(x_i)`, has

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
within that element. Denote the region occupied by element :math:`j` by
:math:`\mathcal D_j` and its volume by
:math:`V_j=\operatorname{vol}(\mathcal D_j)`. If :math:`a_j` is
the element's integrated strength, the spatial PDF is

.. math::

    p_r(\mathbf{r}) = \frac{a_j}{\sum_k a_k}\frac{1}{V_j},
    \qquad \mathbf{r}\in\mathcal D_j.

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
distribution at every radial grid point. Between grid points :math:`u_i` and
:math:`u_{i+1}`, define the local coordinate
:math:`\eta=(u-u_i)/(u_{i+1}-u_i)\in[0,1]`. The conditional energy distribution
is the mixture

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
       +\delta\sin\alpha\cos\alpha\sin\psi
    \end{aligned}

where subscripts :math:`r` and :math:`\alpha` denote partial derivatives. The
coordinate volume element is :math:`R` times this poloidal Jacobian.
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

.. _methods_stellarator_source:

Stellarator Sources
-------------------

Fourier Representation and Volume Jacobian
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The stellarator source uses flux coordinates :math:`(\rho,\theta,\zeta)`,
where :math:`\rho=\sqrt{s}` for normalized toroidal flux :math:`s`,
:math:`\theta` is a poloidal angle, and :math:`\zeta` equals the cylindrical
azimuth. With :math:`N_{\mathrm{fp}}` field periods, define
:math:`\chi_k=m_k\theta-n_kN_{\mathrm{fp}}\zeta`, where :math:`m_k\ge0`
and :math:`n_k` are the integer poloidal and toroidal mode numbers. Their
maximum magnitudes are :math:`m_{\max}=\max_k m_k` and
:math:`|n|_{\max}=\max_k|n_k|`. The surfaces are

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
is exactly :math:`\zeta`.

The geometry must have positive :math:`R` and nested flux surfaces. The sign
of :math:`\tau` depends on the convention for measuring the poloidal angle:
traversing a cross section clockwise or counterclockwise describes the same
surface but gives opposite signs for :math:`\tau`. Physical volume uses its
absolute value, so either convention is allowed. For a valid map, the sign is
consistent throughout the plasma away from the magnetic axis, where the
Jacobian vanishes. OpenMC determines this sign during initialization and sets
:math:`\sigma=+1` for positive :math:`\tau` or :math:`\sigma=-1` for negative
:math:`\tau`. It can then evaluate the nonnegative volume factor as
:math:`J=\sigma R\tau=R|\tau|`. Initialization checks for nonpositive radii
and significant changes in the sign of :math:`\tau` on a setup grid.

The role of :math:`J` is to convert a coordinate element into physical volume:

.. math::

    dV=J(\rho,\theta,\zeta)\,d\rho\,d\theta\,d\zeta.

Multiplying by the emission density :math:`S(\rho)` gives the emission rate
in that element. Define the differential volume and total emission rate as

.. math::

    V'(\rho)=\int_0^{2\pi}\!\int_0^{2\pi}
                J(\rho,\theta,\zeta)\,d\theta\,d\zeta,
    \qquad Q=\int S(u)V'(u)\,du,

where the radial integral covers the tabulated radial range. The normalized
joint PDF with respect to :math:`d\rho\,d\theta\,d\zeta` is therefore

.. math::
    :label: stellarator-source-joint

    p(\rho,\theta,\zeta)
      =\frac{S(\rho)J(\rho,\theta,\zeta)}{Q}.

Integrating out the angles gives the radial marginal. Dividing the joint PDF
by that marginal gives the conditional angular PDF:

.. math::
    :label: stellarator-source-factorization

    p_\rho(\rho)=\frac{S(\rho)V'(\rho)}{Q}, \qquad
    p(\theta,\zeta\mid\rho)
      =\frac{J(\rho,\theta,\zeta)}{V'(\rho)}.

Thus :math:`p(\rho,\theta,\zeta)=p_\rho(\rho)p(\theta,\zeta\mid\rho)`.
OpenMC first samples :math:`\rho` from the radial marginal, then samples the
angles at that radius. The emission profile affects the radial distribution
and cancels from the conditional angular PDF. The following sections derive
the interval polynomials and angular bounds used to sample these two factors.

Interpolation Within a Radial Interval
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

On interval :math:`i`, define its width :math:`h_i=\rho_{i+1}-\rho_i` and local
coordinate :math:`\eta=(\rho-\rho_i)/h_i`. Each Fourier coefficient :math:`C`
and the emission profile are interpolated linearly, using endpoint values
:math:`C_i=C(\rho_i)` and :math:`S_i=S(\rho_i)`:

.. math::

    C(\eta)=(1-\eta)C_i+\eta C_{i+1}, \qquad
    S(\eta)=S_i+\eta(S_{i+1}-S_i).

Functions written with argument :math:`\eta` are evaluated at
:math:`\rho=\rho_i+h_i\eta`. Throughout the following derivation, a single
interval :math:`i` is considered at a time. Polynomial coefficients, control
functions, and cell bounds depend on that interval; the index :math:`i` is
suppressed where it is not written explicitly.

Because Fourier evaluation is linear in the coefficients, at fixed angles
:math:`R`, :math:`Z`, :math:`R_\theta`, and :math:`Z_\theta` are linear in
:math:`\eta`, while

.. math::

    R_\rho=\frac{R_{i+1}-R_i}{h_i}, \qquad
    Z_\rho=\frac{Z_{i+1}-Z_i}{h_i}

are constant within the interval. It follows that :math:`\tau` is linear in
:math:`\eta` and :math:`J=\sigma R\tau` is quadratic. Define the local volume
factor explicitly as

.. math::

    J_i(\eta,\theta,\zeta)=J(\rho_i+h_i\eta,\theta,\zeta).

The differential volume defined above therefore has the interval form

.. math::
    :label: stellarator-source-volume-polynomial

    \begin{aligned}
    V'_i(\eta)
      &\equiv V'(\rho_i+h_i\eta) \\
      &= \int_0^{2\pi}\!\int_0^{2\pi}
         J_i(\eta,\theta,\zeta)\,d\theta\,d\zeta \\
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
these positions, it evaluates :math:`J_i` on a uniform angular grid with sizes

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
      J_i(\eta,\theta_j,\zeta_k).

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

To sample the radial marginal in :eq:`stellarator-source-factorization`,
multiply the linear emission profile by
:eq:`stellarator-source-volume-polynomial`. This gives an unnormalized cubic
radial density with respect to :math:`d\rho`, expressed using the local
coordinate:

.. math::
    :label: stellarator-source-radial-polynomial

    \begin{aligned}
    \widetilde p_i(\eta) &= S(\eta)V'_i(\eta)
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
    M_i &= h_i\int_0^1 \widetilde p_i(\eta)\,d\eta \\
        &= h_i\left(\frac{S_i v_0}{6}
           +\frac{(S_i+S_{i+1})v_m}{3}
           +\frac{S_{i+1}v_1}{6}\right).
    \end{aligned}

The last expression is Simpson's rule, which is exact for a cubic. The masses
sum to the total emission rate :math:`Q=\sum_i M_i`. OpenMC selects interval
:math:`i` with probability :math:`M_i/Q` using an alias table. A common scaling
of all emission values cancels from these probabilities; the implementation
divides them by their maximum before forming the masses and PDFs.

Conditional on choosing interval :math:`i`, the normalized PDF with respect to
:math:`d\eta` is

.. math::

    p(\eta\mid i)=\frac{h_i\widetilde p_i(\eta)}{M_i},
    \qquad 0\le\eta\le1.

The factor :math:`h_i` accounts for the change from :math:`\rho` to
:math:`\eta`; it is constant within the selected interval and cancels from
the rejection-sampling acceptance ratio.

Within the selected interval, OpenMC uses uniform proposals
:math:`\eta=\xi_1` with rejection sampling. An upper bound follows by writing
the cubic in the Bernstein basis:

.. math::

    \widetilde p_i(\eta) = \sum_{j=0}^{3}\beta_j
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
:math:`B_i^{\mathrm{rad}}=\max_j\beta_j` bounds :math:`\widetilde p_i`
throughout the interval.
The implementation adds a small relative margin for floating-point roundoff
and stores the coefficients divided by that bound. A proposal is accepted when

.. math::

    \xi_2 < \frac{\widetilde p_i(\eta)}{B_i^{\mathrm{rad}}}.

On rejection, only :math:`\eta` is redrawn; the selected interval is retained.
Selecting a new interval after a rejection would change its probability in
proportion to its acceptance rate. After acceptance,
:math:`\rho=\rho_i+h_i\eta` has the required radial marginal.

.. _methods_stellarator_angular_sampling:

Conditional Angular Sampling
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

With :math:`\rho` fixed, OpenMC samples the conditional angular PDF in
:eq:`stellarator-source-factorization`. The angular proposal uses a 32 by 32
partition of one field period, with a separate bound for each cell that varies
with radius. To describe the cell geometry, use a toroidal angle
:math:`\zeta_0` restricted to the first field
period and its rescaled counterpart :math:`\varphi`:

.. math::

    \zeta_0\in[0,2\pi/N_{\mathrm{fp}}), \qquad
    \varphi=N_{\mathrm{fp}}\zeta_0\in[0,2\pi).

The geometric volume factor evaluated in these coordinates is

.. math::

    \widetilde J_i(\eta,\theta,\varphi)
      = J_i(\eta,\theta,\zeta_0)
      = J_i(\eta,\theta,\varphi/N_{\mathrm{fp}}).

Both :math:`\theta` and :math:`\varphi` range over :math:`[0,2\pi)`, so each
cell has angular half-width :math:`d=\pi/32` in both coordinates. This cell grid
is separate from the finer quadrature grid used during initialization. The
quadrature grid resolves the density's Fourier coefficients; the cell grid
defines regions for proposing samples. Integrating :math:`\widetilde J_i`
over the rescaled square gives :math:`V'_i(\eta)`: the factor
:math:`d\varphi=N_{\mathrm{fp}}\,d\zeta_0` accounts for all field periods.
The normalized angular PDF in the rescaled square is therefore
:math:`p(\theta,\varphi\mid\rho)=\widetilde J_i/V'_i`.

Within a fixed radial interval, express the quadratic dependence on
:math:`\eta` in the Bernstein basis:

.. math::
    :label: stellarator-source-angular-bernstein

    \begin{aligned}
    \widetilde J_i(\eta,\theta,\varphi)
      &= \sum_{j=0}^2 b_j(\eta)G_j(\theta,\varphi), \\
    b_0(\eta)&=(1-\eta)^2, \qquad
    b_1(\eta)=2\eta(1-\eta), \qquad b_2(\eta)=\eta^2.
    \end{aligned}

Suppressing the interval index on the three angular control functions, their
values follow from the radial evaluations already available during setup:

.. math::

    \begin{aligned}
    G_0 &= \widetilde J_i(0,\theta,\varphi), \\
    G_1 &= 2\widetilde J_i(1/2,\theta,\varphi)
           -\tfrac12[\widetilde J_i(0,\theta,\varphi)+\widetilde J_i(1,\theta,\varphi)], \\
    G_2 &= \widetilde J_i(1,\theta,\varphi).
    \end{aligned}

This is an exact change of polynomial basis. The useful property is that
:math:`b_j\ge0` and :math:`\sum_j b_j=1` on the interval. If a nonnegative
constant :math:`U_{jc}` bounds :math:`G_j` throughout angular cell :math:`c`,
then

.. math::
    :label: stellarator-source-cell-bound

    B_c(\eta)=\sum_{j=0}^2 b_j(\eta)U_{jc}
      \;\ge\; \widetilde J_i(\eta,\theta,\varphi)
      \quad\text{within cell }c.

Thus three bounds computed once per cell produce a bound at every radius in
the interval. For a valid geometry, the control functions :math:`G_j` are
also nonnegative. This follows from the factorization
:math:`\widetilde J_i=R\widetilde\tau`, where
:math:`\widetilde\tau=\sigma\tau` is the nonnegative poloidal Jacobian.
Both factors are linear in :math:`\eta`. Writing their endpoint values as
:math:`R_i,R_{i+1}` and :math:`\widetilde\tau_i,\widetilde\tau_{i+1}` gives

.. math::

    \begin{aligned}
    G_0 &= R_i\widetilde\tau_i, \\
    G_1 &= \tfrac12(R_i\widetilde\tau_{i+1}+R_{i+1}\widetilde\tau_i), \\
    G_2 &= R_{i+1}\widetilde\tau_{i+1}.
    \end{aligned}

Each expression is nonnegative because the endpoint radii are positive and
the endpoint values of :math:`\widetilde\tau` are nonnegative. The
:math:`G_j` are unnormalized angular functions; the cell-selection
probabilities are constructed from their bounds below.

Bounding the Density on Angular Cells
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Each :math:`G_j` is a finite Fourier series in the two angular coordinates:

.. math::

    G_j(\theta,\varphi)
      = \sum_{m=-M}^{M}\sum_{n=-L}^{L}
        \widehat G_{j,mn}e^{\mathrm{i}(m\theta+n\varphi)},
    \qquad M=3m_{\max},\quad L=3|n|_{\max}.

Initialization recovers these coefficients using separable discrete Fourier
sums on the quadrature grid. Resolving all coefficients requires the stronger
conditions :math:`N_\theta>2M` and :math:`N_\zeta>2L`, which the grid sizes
above also satisfy. The implementation evaluates the sums directly and needs
no external FFT library.

Let :math:`\boldsymbol{a}_c=(\theta_c,\varphi_c)` be a cell center. Any point
in the cell can be written as :math:`\boldsymbol{a}_c+\boldsymbol{v}`, where
:math:`|v_\theta|,|v_\varphi|\le d`. The goal is to bound the largest possible
increase from :math:`G_j(\boldsymbol{a}_c)` over all such displacements. Taylor's
theorem separates this increase into a contribution from the slope at the
center and a remainder accounting for curvature.

For a fixed displacement, reduce the problem to one dimension by following
the straight line from the center to the point:

.. math::

    g(\lambda)=G_j(\boldsymbol{a}_c+\lambda\boldsymbol{v}),
    \qquad 0\le\lambda\le1.

Here :math:`\lambda` measures progress along the line; the radial coordinate
:math:`\eta` remains fixed during angular sampling. Taylor's theorem with the
Lagrange remainder gives

.. math::
    :label: stellarator-source-line-taylor

    g(1)=g(0)+g'(0)+\frac12 g''(\lambda_*),
    \qquad \lambda_*\in(0,1).

This is an exact identity. The second derivative is evaluated at an unknown
intermediate point, and its term accounts for the entire departure from the
linear approximation. Bounding this derivative along the line therefore bounds
the full remainder, including effects that a longer Taylor expansion would
express through higher derivatives. No additional terms need to be neglected.

The center value and slope follow from the chain rule:

.. math::

    \begin{aligned}
    g(0)&=G_j(\boldsymbol{a}_c), \\
    g'(0)&=v_\theta\partial_\theta G_j(\boldsymbol{a}_c)
            +v_\varphi\partial_\varphi G_j(\boldsymbol{a}_c).
    \end{aligned}

Using the maximum displacement in each coordinate bounds the slope term by

.. math::

    g'(0)\le d\left(|\partial_\theta G_j(\boldsymbol{a}_c)|
                    +|\partial_\varphi G_j(\boldsymbol{a}_c)|\right).

To bound the remainder, differentiate the Fourier series twice along the line:

.. math::

    \begin{aligned}
    g''(\lambda)
      &= -\sum_{m,n}\widehat G_{j,mn}(mv_\theta+nv_\varphi)^2 \\
      &\qquad\times e^{\mathrm{i}[m(\theta_c+\lambda v_\theta)
                                  +n(\varphi_c+\lambda v_\varphi)]}.
    \end{aligned}

Each exponential has magnitude one. The triangle inequality then gives a
bound valid at every point along the line:

.. math::

    |g''(\lambda)|\le\sum_{m,n}|\widehat G_{j,mn}|
                                      |mv_\theta+nv_\varphi|^2.

Because the point lies within the cell,

.. math::

    |mv_\theta+nv_\varphi|
      \le |m|\,|v_\theta|+|n|\,|v_\varphi|
      \le d(|m|+|n|).

Consequently,

.. math::

    |g''(\lambda)|\le d^2 C_j, \qquad
    C_j=\sum_{m,n}|\widehat G_{j,mn}|(|m|+|n|)^2.

The quantity :math:`C_j` depends only on the Fourier coefficients and bounds
curvature for every cell and displacement. Substituting this result and the
slope bound into :eq:`stellarator-source-line-taylor` gives

.. math::
    :label: stellarator-source-cell-taylor

    \begin{aligned}
    G_j(\boldsymbol{a}_c+\boldsymbol{v})
    &\le G_j(\boldsymbol{a}_c)
       +d\left(|\partial_\theta G_j(\boldsymbol{a}_c)|
              +|\partial_\varphi G_j(\boldsymbol{a}_c)|\right) \\
    &\qquad +\frac{d^2}{2}C_j.
    \end{aligned}

Every quantity on the right-hand side is available during initialization and
is independent of the unknown displacement. The inequality therefore bounds
the function throughout the cell, including peaks between evaluation points.
The curvature allowance remains present even when the center gradient is zero.
The center value and its first derivatives retain the cancellations among
Fourier modes, allowing a tighter bound than separately bounding all factors of
the Jacobian. Taking absolute values in the slope and remainder terms makes
the bound conservative; it need not equal a value attained inside the cell.

For each control function :math:`G_j` and cell :math:`c`, OpenMC evaluates the
right-hand side of :eq:`stellarator-source-cell-taylor` using the function's
value and slopes at the cell center and the curvature bound :math:`C_j`.
All Fourier coefficients contribute to these calculations, so the bound
accounts for variation inside the cell as well as at its center. The resulting
bounds :math:`U_{jc}` are stored during initialization and combined at the
sampled radius as described below.

Selecting and Sampling Angular Cells
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For each control function, define the sum of its cell bounds

.. math::

    T_j=\sum_c U_{jc}.

Initialization builds three alias tables, one for each distribution
:math:`U_{jc}/T_j`. At the sampled radius, OpenMC selects a control index
:math:`j` with probability

.. math::

    P(j\mid\eta)=\frac{b_j(\eta)T_j}{\sum_k b_k(\eta)T_k},

then selects a cell from that control's alias table. Summing over control
indices shows that the resulting cell probability is

.. math::

    P(c\mid\eta)
      =\sum_{j:T_j>0} P(j\mid\eta)\frac{U_{jc}}{T_j}
      =\frac{B_c(\eta)}{\sum_{c'}B_{c'}(\eta)}.

Components with :math:`T_j=0` have zero selection probability and require no
alias table. This mixture samples the radius-dependent bounds using three
fixed tables, avoiding construction of a new table for each sampled radius.

The sampler proposes a point uniformly inside the chosen cell and accepts it
when

.. math::

    \xi B_c(\eta) < \widetilde J_i(\eta,\theta,\varphi).

All cells have the same area, so the proposal density within a cell is
proportional to :math:`B_c`. Multiplying it by the acceptance probability
:math:`\widetilde J_i/B_c` gives a density proportional to
:math:`\widetilde J_i`, as required.
On rejection, both the control index and the cell are redrawn while retaining
the radial interval and :math:`\rho`.

After acceptance, an independently and uniformly chosen field-period index
:math:`\ell\in\{0,\ldots,N_{\mathrm{fp}}-1\}` gives the full-torus angle

.. math::

    \zeta=\zeta_0+\frac{2\pi\ell}{N_{\mathrm{fp}}}
         =\frac{\varphi+2\pi\ell}{N_{\mathrm{fp}}}.

The sampled :math:`R` and :math:`Z` apply in every field period. The cell
partition affects rejection efficiency while the accepted points follow the
continuous Fourier density. The Taylor inequality supplies bounds in exact
arithmetic; a runtime check terminates sampling if a density exceeds its cell
bound. Nonfinite or negative densities and excessive rejection also produce
errors.

Sampling Sequence and Approximation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The steps in the stellarator source's sampling method are:

1. Select a radial interval using its integrated mass :math:`M_i`.
2. Sample the local radial coordinate by cubic-PDF rejection within that
   interval, and calculate :math:`\rho`.
3. At that fixed radius, select a control index and an angular cell, propose
   angles uniformly inside the cell, and reject using the local bound. Repeat
   angular proposals until accepted, then choose a field period uniformly.
4. Map the accepted position to :math:`(R\cos\zeta,R\sin\zeta,Z)`.
5. Sample an isotropic flight direction.
6. Sample energy from the common spectrum or the mixture of spectra bracketing
   :math:`\rho`, then sample emission time and apply their weight factors.

Radial integration and rejection sampling are consistent with the piecewise
linear coefficient and emission tables, up to floating-point effects. Those
tables approximate the underlying equilibrium and emission profile; increasing
the number of radial surfaces can reduce interpolation error. The finite
Fourier representation also limits the geometric detail retained. The
Bernstein representation and angular-cell partition introduce no additional
approximation to this density; the angular bounds control rejection efficiency.

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

    p_{\mathcal{A}}(\boldsymbol{x}) =
    \frac{p(\boldsymbol{x})\,\mathbf{1}_{\mathcal{A}}(\boldsymbol{x})}
         {\int_{\mathcal{A}}p(\boldsymbol{y})\,d\boldsymbol{y}},

where :math:`\boldsymbol{x}=(\mathbf{r},\mathbf{\Omega},E,t)` denotes a
phase-space point and :math:`\mathbf{1}_{\mathcal A}` is one inside the
acceptance region and zero outside it. The base source wrapper repeats
sampling until a candidate passes, or, with the ``kill`` strategy, returns a
failing candidate with weight zero. Independent sources handle constraints
internally: they resample positions for spatial constraints and energies for
energy constraints and the available energy range. Their time samples are drawn
directly. Constraints on an element source in a
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
:math:`w\leftarrow w\,p(\boldsymbol{x})/g(\boldsymbol{x})`. The proposal must
cover every region where the desired distribution contributes. When multiple
independent factors are biased, their corrections multiply. Mesh and point-cloud
distributions can similarly change their discrete selection probabilities and
apply the ratio of the original to the biased normalized probability. The
derivation and available biasing mechanisms are described in
:ref:`methods_source_biasing`.

.. [TrefethenWeideman] Lloyd N. Trefethen and J. A. C. Weideman,
   `The Exponentially Convergent Trapezoidal Rule
   <https://people.maths.ox.ac.uk/trefethen/sirev56-3_385.pdf>`_,
   *SIAM Review*, **56** (3), 385--458 (2014), Corollary 3.3, p. 398.
