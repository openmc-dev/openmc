.. _methods_cmfd:

================================================================
Nonlinear Diffusion Acceleration - Coarse Mesh Finite Difference
================================================================

This page section discusses how nonlinear diffusion acceleration (NDA) using coarse mesh finite difference (CMFD) is implemented into OpenMC. Before we get into the theory, general notation for this section is discussed.

--------
Notation
--------

Before deriving NDA relationships, notation is explained. If a parameter has a
:math:`\overline{\cdot}`, it is surface area-averaged and if it has a
:math:`\overline{\overline\cdot}`, it is volume-averaged. When describing a
specific cell in the geometry, indices :math:`(i,j,k)` are used which correspond
to directions :math:`(x,y,z)`. In most cases, the same operation is performed in
all three directions. To compactly write this, an arbitrary direction set
:math:`(u,v,w)` that corresponds to cell indices :math:`(l,m,n)` is used. Note
that :math:`u` and :math:`l` do not have to correspond to :math:`x` and
:math:`i`. However, if :math:`u` and :math:`l` correspond to :math:`y` and
:math:`j`, :math:`v` and :math:`w` correspond to :math:`x` and :math:`z`
directions. An example of this is shown in the following expression:

.. math::
    :label: not1

    \sum\limits_{u\in(x,y,z)}\left\langle\overline{J}^{u,g}_{l+1/2,m,n}
    \Delta_m^v\Delta_n^w\right\rangle 

Here, :math:`u` takes on each direction one at a time. The parameter :math:`J`
is surface area-averaged over the transverse indices :math:`m` and :math:`n`
located at :math:`l+1/2`.  Usually, spatial indices are listed as subscripts and
the direction as a superscript. Energy group indices represented by :math:`g`
and :math:`h` are also listed as superscripts here. The group :math:`g` is the
group of interest and, if present, :math:`h` is all groups. Finally, any
parameter surrounded by :math:`\left\langle\cdot\right\rangle` represents a
tally quantity that can be edited from an MC solution.

------
Theory
------

NDA is a diffusion model that has equivalent physics to a transport model. There
are many different methods that can be classified as NDA. The CMFD method is a
type of NDA that represents second order multigroup diffusion equations on a
coarse spatial mesh.  Whether a transport model or diffusion model is used to
represent the distribution of neutrons, these models must satisfy the *neutron
balance equation*. This balance is represented by the following formula for a
specific energy group :math:`g` in cell :math:`(l,m,n)`:

.. math::
    :label: eq_neut_bal

    \sum\limits_{u\in(x,y,z)}\left(\left\langle\overline{J}^{u,g}_{l+1/2,m,n}
    \Delta_m^v\Delta_n^w\right\rangle -
    \left\langle\overline{J}^{u,g}_{l-1/2,m,n}
    \Delta_m^v\Delta_n^w\right\rangle\right)
    +
    \left\langle\overline{\overline\Sigma}_{t_{l,m,n}}^g
    \overline{\overline\phi}_{l,m,n}^g\Delta_l^u\Delta_m^v\Delta_n^w\right\rangle
    = \\
    \sum\limits_{h=1}^G\left\langle
    \overline{\overline{\nu_s\Sigma}}_{s_{l,m,n}}^{h\rightarrow
    g}\overline{\overline\phi}_{l,m,n}^h\Delta_l^u\Delta_m^v\Delta_n^w
    \right\rangle
    +
    \frac{1}{k_{eff}}\sum\limits_{h=1}^G
    \left\langle\overline{\overline{\nu_f\Sigma}}_{f_{l,m,n}}^{h\rightarrow
    g}\overline{\overline\phi}_{l,m,n}^h
    \Delta_l^u\Delta_m^v\Delta_n^w\right\rangle.

In eq. :eq:`eq_neut_bal` the parameters are defined as:

* :math:`\left\langle\overline{J}^{u,g}_{l\pm
  1/2,m,n}\Delta_m^v\Delta_n^w\right\rangle` --- surface area-integrated net
  current over surface :math:`(l\pm 1/2,m,n)` with surface normal in direction
  $u$ in energy group :math:`g`. By dividing this quantity by the transverse
  area, :math:`\Delta_m^v\Delta_n^w`, the surface area-averaged net current can
  be computed.
* :math:`\left\langle\overline{\overline\Sigma}_{t_{l,m,n}}^g
  \overline{\overline\phi}_{l,m,n}^g\Delta_l^u\Delta_m^v\Delta_n^w\right\rangle`
  --- volume-integrated total reaction rate over energy group :math:`g`.
  * :math:`\left\langle\overline{\overline{\nu_s\Sigma}}_{s_{l,m,n}}^{h\rightarrow
  g}
  \overline{\overline\phi}_{l,m,n}^h\Delta_l^u\Delta_m^v\Delta_n^w\right\rangle`
  --- volume-integrated scattering production rate of neutrons that begin with
  energy in group :math:`h` and exit reaction in group :math:`g`. This reaction
  rate also includes the energy transfer of reactions (except fission) that
  produce multiple neutrons such as (n, 2n); hence, the need for :math:`\nu_s`
  to represent neutron multiplicity.
* :math:`k_{eff}` --- core multiplication factor.
* :math:`\left\langle\overline{\overline{\nu_f\Sigma}}_{f_{l,m,n}}^{h\rightarrow
  g}\overline{\overline\phi}_{l,m,n}^h\Delta_l^u\Delta_m^v\Delta_n^w\right\rangle`
  --- volume-integrated fission production rate of neutrons from fissions in
  group :math:`h` that exit in group :math:`g`.

Each quantity in :math:`\left\langle\cdot\right\rangle` represents a scalar value that
is obtained from an MC tally. A good verification step when using an MC is
to make sure that tallies satisfy this balance equation within statistics. No
NDA acceleration can be performed if the balance equation is not satisfied.

There are three major steps to consider when performing NDA: (1) calculation of
macroscopic cross sections and nonlinear parameters, (2) solving an eigenvalue
problem with a system of linear equations, and (3) modifying MC source
distribution to align with the NDA solution on a chosen mesh. This process is
illustrated as a flow chart below. After a batch of neutrons
is simulated, NDA can take place. Each of the steps described above is described
in detail in the following sections.

.. tikz:: Flow chart of NDA process 
   :libs: shapes, snakes, shadows, arrows, calc, decorations.markings, patterns, fit, matrix, spy
   :include: cmfd_tikz/cmfd_flow.tikz
