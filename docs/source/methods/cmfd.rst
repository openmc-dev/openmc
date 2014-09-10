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
