.. _io_tallies:

====================================
Tallies Specification -- tallies.xml
====================================

The tallies.xml file allows the user to tell the code what results he/she is
interested in, e.g. the fission rate in a given cell or the current across a
given surface. There are two pieces of information that determine what
quantities should be scored. First, one needs to specify what region of phase
space should count towards the tally and secondly, the actual quantity to be
scored also needs to be specified. The first set of parameters we call *filters*
since they effectively serve to filter events, allowing some to score and
preventing others from scoring to the tally.

The structure of tallies in OpenMC is flexible in that any combination of
filters can be used for a tally. The following types of filter are available:
cell, universe, material, surface, birth region, pre-collision energy,
post-collision energy, and an arbitrary structured mesh.

The three valid elements in the tallies.xml file are ``<tally>``, ``<mesh>``,
and ``<assume_separate>``.

.. _tally:

-------------------
``<tally>`` Element
-------------------

The ``<tally>`` element accepts the following sub-elements:

  :name:
    An optional string name to identify the tally in summary output
    files. This string is limited to 52 characters for formatting purposes.

    *Default*: ""

  :filter:
    Specify a filter that modifies tally behavior. Most tallies (e.g. ``cell``,
    ``energy``, and ``material``) restrict the tally so that only particles
    within certain regions of phase space contribute to the tally.  Others
    (e.g. ``delayedgroup`` and ``energyfunction``) can apply some other function
    to the scored values. This element and its attributes/sub-elements are
    described below.

    .. note::
        You may specify zero, one, or multiple filters to apply to the tally. To
        specify multiple filters, you must use multiple ``<filter>`` elements.

    The ``filter`` element has the following attributes/sub-elements:

      :type:
        The type of the filter. Accepted options are "cell", "cellborn",
        "material", "universe", "energy", "energyout", "mu", "polar",
        "azimuthal", "mesh", "distribcell", "delayedgroup", and
        "energyfunction".

      :bins:
         A description of the bins for each type of filter can be found in
         :ref:`filter_types`.

      :energy:
        ``energyfunction`` filters multiply tally scores by an arbitrary
        function. The function is described by a piecewise linear-linear set of
        (energy, y) values. This entry specifies the energy values. The function
        will be evaluated as zero outside of the bounds of this energy grid.
        (Only used for ``energyfunction`` filters)

      :y:
        ``energyfunction`` filters multiply tally scores by an arbitrary
        function. The function is described by a piecewise linear-linear set of
        (energy, y) values. This entry specifies the y values. (Only used
        for ``energyfunction`` filters)

  :nuclides:
    If specified, the scores listed will be for particular nuclides, not the
    summation of reactions from all nuclides. The format for nuclides should be
    [Atomic symbol]-[Mass number], e.g. "U-235". The reaction rate for all
    nuclides can be obtained with "total". For example, to obtain the reaction
    rates for U-235, Pu-239, and all nuclides in a material, this element should
    be:

    .. code-block:: xml

        <nuclides>U-235 Pu-239 total</nuclides>

    *Default*: total

  :estimator:
    The estimator element is used to force the use of either ``analog``,
    ``collision``, or ``tracklength`` tally estimation.  ``analog`` is generally
    the least efficient though it can be used with every score type.
    ``tracklength`` is generally the most efficient, but neither ``tracklength``
    nor ``collision`` can be used to score a tally that requires post-collision
    information.  For example, a scattering tally with outgoing energy filters
    cannot be used with ``tracklength`` or ``collision`` because the code will
    not know the outgoing energy distribution.

    *Default*: ``tracklength`` but will revert to ``analog`` if necessary.

  :scores:
    A space-separated list of the desired responses to be accumulated. The accepted
    options are listed in the following tables:

    .. table:: **Flux scores: units are particle-cm per source particle.**

        +----------------------+---------------------------------------------------+
        |Score                 | Description                                       |
        +======================+===================================================+
        |flux                  |Total flux.                                        |
        +----------------------+---------------------------------------------------+
        |flux-YN               |Spherical harmonic expansion of the direction of   |
        |                      |motion :math:`\left(\Omega\right)` of the total    |
        |                      |flux. This score will tally all of the harmonic    |
        |                      |moments of order 0 to N.  N must be between 0 and  |
        |                      |10.                                                |
        +----------------------+---------------------------------------------------+

    .. table:: **Reaction scores: units are reactions per source particle.**

        +----------------------+---------------------------------------------------+
        |Score                 | Description                                       |
        +======================+===================================================+
        |absorption            |Total absorption rate. This accounts for all       |
        |                      |reactions which do not produce secondary neutrons  |
        |                      |as well as fission.                                |
        +----------------------+---------------------------------------------------+
        |elastic               |Elastic scattering reaction rate.                  |
        +----------------------+---------------------------------------------------+
        |fission               |Total fission reaction rate.                       |
        +----------------------+---------------------------------------------------+
        |scatter               |Total scattering rate. Can also be identified with |
        |                      |the "scatter-0" response type.                     |
        +----------------------+---------------------------------------------------+
        |scatter-N             |Tally the N\ :sup:`th` \ scattering moment, where N|
        |                      |is the Legendre expansion order of the change in   |
        |                      |particle angle :math:`\left(\mu\right)`. N must be |
        |                      |between 0 and 10. As an example, tallying the 2\   |
        |                      |:sup:`nd` \ scattering moment would be specified as|
        |                      |``<scores>scatter-2</scores>``.                    |
        +----------------------+---------------------------------------------------+
        |scatter-PN            |Tally all of the scattering moments from order 0 to|
        |                      |N, where N is the Legendre expansion order of the  |
        |                      |change in particle angle                           |
        |                      |:math:`\left(\mu\right)`. That is, "scatter-P1" is |
        |                      |equivalent to requesting tallies of "scatter-0" and|
        |                      |"scatter-1".  Like for "scatter-N", N must be      |
        |                      |between 0 and 10. As an example, tallying up to the|
        |                      |2\ :sup:`nd` \ scattering moment would be specified|
        |                      |as ``<scores> scatter-P2 </scores>``.              |
        +----------------------+---------------------------------------------------+
        |scatter-YN            |"scatter-YN" is similar to "scatter-PN" except an  |
        |                      |additional expansion is performed for the incoming |
        |                      |particle direction :math:`\left(\Omega\right)`     |
        |                      |using the real spherical harmonics.  This is useful|
        |                      |for performing angular flux moment weighting of the|
        |                      |scattering moments. Like "scatter-PN", "scatter-YN"|
        |                      |will tally all of the moments from order 0 to N; N |
        |                      |again must be between 0 and 10.                    |
        +----------------------+---------------------------------------------------+
        |total                 |Total reaction rate.                               |
        +----------------------+---------------------------------------------------+
        |total-YN              |The total reaction rate expanded via spherical     |
        |                      |harmonics about the direction of motion of the     |
        |                      |neutron, :math:`\Omega`. This score will tally all |
        |                      |of the harmonic moments of order 0 to N.  N must be|
        |                      |between 0 and 10.                                  |
        +----------------------+---------------------------------------------------+
        |(n,2nd)               |(n,2nd) reaction rate.                             |
        +----------------------+---------------------------------------------------+
        |(n,2n)                |(n,2n) reaction rate.                              |
        +----------------------+---------------------------------------------------+
        |(n,3n)                |(n,3n) reaction rate.                              |
        +----------------------+---------------------------------------------------+
        |(n,na)                |(n,n\ :math:`\alpha`\ ) reaction rate.             |
        +----------------------+---------------------------------------------------+
        |(n,n3a)               |(n,n3\ :math:`\alpha`\ ) reaction rate.            |
        +----------------------+---------------------------------------------------+
        |(n,2na)               |(n,2n\ :math:`\alpha`\ ) reaction rate.            |
        +----------------------+---------------------------------------------------+
        |(n,3na)               |(n,3n\ :math:`\alpha`\ ) reaction rate.            |
        +----------------------+---------------------------------------------------+
        |(n,np)                |(n,np) reaction rate.                              |
        +----------------------+---------------------------------------------------+
        |(n,n2a)               |(n,n2\ :math:`\alpha`\ ) reaction rate.            |
        +----------------------+---------------------------------------------------+
        |(n,2n2a)              |(n,2n2\ :math:`\alpha`\ ) reaction rate.           |
        +----------------------+---------------------------------------------------+
        |(n,nd)                |(n,nd) reaction rate.                              |
        +----------------------+---------------------------------------------------+
        |(n,nt)                |(n,nt) reaction rate.                              |
        +----------------------+---------------------------------------------------+
        |(n,nHe-3)             |(n,n\ :sup:`3`\ He) reaction rate.                 |
        +----------------------+---------------------------------------------------+
        |(n,nd2a)              |(n,nd2\ :math:`\alpha`\ ) reaction rate.           |
        +----------------------+---------------------------------------------------+
        |(n,nt2a)              |(n,nt2\ :math:`\alpha`\ ) reaction rate.           |
        +----------------------+---------------------------------------------------+
        |(n,4n)                |(n,4n) reaction rate.                              |
        +----------------------+---------------------------------------------------+
        |(n,2np)               |(n,2np) reaction rate.                             |
        +----------------------+---------------------------------------------------+
        |(n,3np)               |(n,3np) reaction rate.                             |
        +----------------------+---------------------------------------------------+
        |(n,n2p)               |(n,n2p) reaction rate.                             |
        +----------------------+---------------------------------------------------+
        |(n,n*X*)              |Level inelastic scattering reaction rate. The *X*  |
        |                      |indicates what which inelastic level, e.g., (n,n3) |
        |                      |is third-level inelastic scattering.               |
        +----------------------+---------------------------------------------------+
        |(n,nc)                |Continuum level inelastic scattering reaction rate.|
        +----------------------+---------------------------------------------------+
        |(n,gamma)             |Radiative capture reaction rate.                   |
        +----------------------+---------------------------------------------------+
        |(n,p)                 |(n,p) reaction rate.                               |
        +----------------------+---------------------------------------------------+
        |(n,d)                 |(n,d) reaction rate.                               |
        +----------------------+---------------------------------------------------+
        |(n,t)                 |(n,t) reaction rate.                               |
        +----------------------+---------------------------------------------------+
        |(n,3He)               |(n,\ :sup:`3`\ He) reaction rate.                  |
        +----------------------+---------------------------------------------------+
        |(n,a)                 |(n,\ :math:`\alpha`\ ) reaction rate.              |
        +----------------------+---------------------------------------------------+
        |(n,2a)                |(n,2\ :math:`\alpha`\ ) reaction rate.             |
        +----------------------+---------------------------------------------------+
        |(n,3a)                |(n,3\ :math:`\alpha`\ ) reaction rate.             |
        +----------------------+---------------------------------------------------+
        |(n,2p)                |(n,2p) reaction rate.                              |
        +----------------------+---------------------------------------------------+
        |(n,pa)                |(n,p\ :math:`\alpha`\ ) reaction rate.             |
        +----------------------+---------------------------------------------------+
        |(n,t2a)               |(n,t2\ :math:`\alpha`\ ) reaction rate.            |
        +----------------------+---------------------------------------------------+
        |(n,d2a)               |(n,d2\ :math:`\alpha`\ ) reaction rate.            |
        +----------------------+---------------------------------------------------+
        |(n,pd)                |(n,pd) reaction rate.                              |
        +----------------------+---------------------------------------------------+
        |(n,pt)                |(n,pt) reaction rate.                              |
        +----------------------+---------------------------------------------------+
        |(n,da)                |(n,d\ :math:`\alpha`\ ) reaction rate.             |
        +----------------------+---------------------------------------------------+
        |*Arbitrary integer*   |An arbitrary integer is interpreted to mean the    |
        |                      |reaction rate for a reaction with a given ENDF MT  |
        |                      |number.                                            |
        +----------------------+---------------------------------------------------+

    .. table:: **Particle production scores: units are particles produced per
               source particles.**

        +----------------------+---------------------------------------------------+
        |Score                 | Description                                       |
        +======================+===================================================+
        |delayed-nu-fission    |Total production of delayed neutrons due to        |
        |                      |fission.                                           |
        +----------------------+---------------------------------------------------+
        |prompt-nu-fission     |Total production of prompt neutrons due to         |
        |                      |fission.                                           |
        +----------------------+---------------------------------------------------+
        |nu-fission            |Total production of neutrons due to fission.       |
        +----------------------+---------------------------------------------------+
        |nu-scatter,           |These scores are similar in functionality to their |
        |nu-scatter-N,         |``scatter*`` equivalents except the total          |
        |nu-scatter-PN,        |production of neutrons due to scattering is scored |
        |nu-scatter-YN         |vice simply the scattering rate. This accounts for |
        |                      |multiplicity from (n,2n), (n,3n), and (n,4n)       |
        |                      |reactions.                                         |
        +----------------------+---------------------------------------------------+

    .. table:: **Miscellaneous scores: units are indicated for each.**

        +----------------------+---------------------------------------------------+
        |Score                 | Description                                       |
        +======================+===================================================+
        |current               |Partial currents on the boundaries of each cell in |
        |                      |a mesh. Units are particles per source             |
        |                      |particle. Note that this score can only be used if |
        |                      |a mesh filter has been specified. Furthermore, it  |
        |                      |may not be used in conjunction with any other      |
        |                      |score.                                             |
        +----------------------+---------------------------------------------------+
        |events                |Number of scoring events. Units are events per     |
        |                      |source particle.                                   |
        +----------------------+---------------------------------------------------+
        |inverse-velocity      |The flux-weighted inverse velocity where the       |
        |                      |velocity is in units of centimeters per second.    |
        +----------------------+---------------------------------------------------+
        |kappa-fission         |The recoverable energy production rate due to      |
        |                      |fission. The recoverable energy is defined as the  |
        |                      |fission product kinetic energy, prompt and delayed |
        |                      |neutron kinetic energies, prompt and delayed       |
        |                      |:math:`\gamma`-ray total energies, and the total   |
        |                      |energy released by the delayed :math:`\beta`       |
        |                      |particles. The neutrino energy does not contribute |
        |                      |to this response. The prompt and delayed           |
        |                      |:math:`\gamma`-rays are assumed to deposit their   |
        |                      |energy locally. Units are eV per source particle.  |
        +----------------------+---------------------------------------------------+
        |fission-q-prompt      |The prompt fission energy production rate. This    |
        |                      |energy comes in the form of fission fragment       |
        |                      |nuclei, prompt neutrons, and prompt                |
        |                      |:math:`\gamma`-rays. This value depends on the     |
        |                      |incident energy and it requires that the nuclear   |
        |                      |data library contains the optional fission energy  |
        |                      |release data. Energy is assumed to be deposited    |
        |                      |locally. Units are eV per source particle.         |
        +----------------------+---------------------------------------------------+
        |fission-q-recoverable |The recoverable fission energy production rate.    |
        |                      |This energy comes in the form of fission fragment  |
        |                      |nuclei, prompt and delayed neutrons, prompt and    |
        |                      |delayed :math:`\gamma`-rays, and delayed           |
        |                      |:math:`\beta`-rays. This tally differs from the    |
        |                      |kappa-fission tally in that it is dependent on     |
        |                      |incident neutron energy and it requires that the   |
        |                      |nuclear data library contains the optional fission |
        |                      |energy release data. Energy is assumed to be       |
        |                      |deposited locally. Units are eV per source         |
        |                      |paticle.                                           |
        +----------------------+---------------------------------------------------+
        |decay-rate            |The delayed-nu-fission-weighted decay rate where   |
        |                      |the decay rate is in units of inverse seconds.     |
        +----------------------+---------------------------------------------------+

    .. note::
       The ``analog`` estimator is actually identical to the ``collision``
       estimator for the flux and inverse-velocity scores.

  :trigger:
    Precision trigger applied to all filter bins and nuclides for this tally.
    It must specify the trigger's type, threshold and scores to which it will
    be applied. It has the following attributes/sub-elements:

   :type:
     The type of the trigger. Accepted options are "variance", "std_dev",
     and "rel_err".

     :variance:
       Variance of the batch mean :math:`\sigma^2`

     :std_dev:
       Standard deviation of the batch mean :math:`\sigma`

     :rel_err:
       Relative error of the batch mean :math:`\frac{\sigma}{\mu}`

     *Default*: None

   :threshold:
     The precision trigger's convergence criterion for tallied values.

     *Default*: None

   :scores:
     The score(s) in this tally to which the trigger should be applied.

     .. note:: The ``scores`` in ``trigger`` must have been defined in
               ``scores`` in ``tally``. An optional "all" may be used to
               select all scores in this tally.

     *Default*: "all"

  :derivative:
    The id of a ``derivative`` element. This derivative will be applied to all
    scores in the tally. Differential tallies are currently only implemented
    for collision and analog estimators.

     *Default*: None

.. _filter_types:

Filter Types
++++++++++++

For each filter type, the following table describes what the ``bins`` attribute
should be set to:

:cell:
  A list of unique IDs for cells in which the tally should be accumulated.

:cellborn:
  This filter allows the tally to be scored to only when particles were
  originally born in a specified cell. A list of cell IDs should be given.

:material:
  A list of unique IDs for matreials in which the tally should be accumulated.

:universe:
  A list of unique IDs for universes in which the tally should be accumulated.

:energy:
  In continuous-energy mode, this filter should be provided as a
  monotonically increasing list of bounding **pre-collision** energies
  for a number of groups. For example, if this filter is specified as

  .. code-block:: xml

      <filter type="energy" bins="0.0 1.0e6 20.0e6" />

  then two energy bins will be created, one with energies between 0 and
  1 MeV and the other with energies between 1 and 20 MeV.

  In multi-group mode the bins provided must match group edges
  defined in the multi-group library.

:energyout:
  In continuous-energy mode, this filter should be provided as a
  monotonically increasing list of bounding **post-collision** energies
  for a number of groups. For example, if this filter is specified as

  .. code-block:: xml

      <filter type="energyout" bins="0.0 1.0e6 20.0e6" />

  then two post-collision energy bins will be created, one with
  energies between 0 and 1 MeV and the other with energies between
  1 and 20 MeV.

  In multi-group mode the bins provided must match group edges
  defined in the multi-group library.

:mu:
  A monotonically increasing list of bounding **post-collision** cosines
  of the change in a particle's angle (i.e., :math:`\mu = \hat{\Omega}
  \cdot \hat{\Omega}'`), which represents a portion of the possible
  values of :math:`[-1,1]`.  For example, spanning all of :math:`[-1,1]`
  with five equi-width bins can be specified as:

  .. code-block:: xml

      <filter type="mu" bins="-1.0 -0.6 -0.2 0.2 0.6 1.0" />

  Alternatively, if only one value is provided as a bin, OpenMC will
  interpret this to mean the complete range of :math:`[-1,1]` should
  be automatically subdivided in to the provided value for the bin.
  That is, the above example of five equi-width bins spanning
  :math:`[-1,1]` can be instead written as:

  .. code-block:: xml

      <filter type="mu" bins="5" />

:polar:
  A monotonically increasing list of bounding particle polar angles
  which represents a portion of the possible values of :math:`[0,\pi]`.
  For example, spanning all of :math:`[0,\pi]` with five equi-width
  bins can be specified as:

  .. code-block:: xml

      <filter type="polar" bins="0.0 0.6283 1.2566 1.8850 2.5132 3.1416"/>

  Alternatively, if only one value is provided as a bin, OpenMC will
  interpret this to mean the complete range of :math:`[0,\pi]` should
  be automatically subdivided in to the provided value for the bin.
  That is, the above example of five equi-width bins spanning
  :math:`[0,\pi]` can be instead written as:

  .. code-block:: xml

      <filter type="polar" bins="5" />

:azimuthal:
  A monotonically increasing list of bounding particle azimuthal angles
  which represents a portion of the possible values of :math:`[-\pi,\pi)`.
  For example, spanning all of :math:`[-\pi,\pi)` with two equi-width
  bins can be specified as:

  .. code-block:: xml

      <filter type="azimuthal" bins="0.0 3.1416 6.2832" />

  Alternatively, if only one value is provided as a bin, OpenMC will
  interpret this to mean the complete range of :math:`[-\pi,\pi)` should
  be automatically subdivided in to the provided value for the bin.
  That is, the above example of five equi-width bins spanning
  :math:`[-\pi,\pi)` can be instead written as:

  .. code-block:: xml

      <filter type="azimuthal" bins="2" />

:mesh:
  The unique ID of a structured mesh to be tallied over.

:distribcell:
  The single cell which should be tallied uniquely for all instances.

  .. note:: The distribcell filter will take a single cell ID and will tally
            each unique occurrence of that cell separately. This filter will not
            accept more than one cell ID. It is not recommended to combine this
            filter with a cell or mesh filter.

:delayedgroup:
  A list of delayed neutron precursor groups for which the tally should
  be accumulated. For instance, to tally to all 6 delayed groups in the
  ENDF/B-VII.1 library the filter is specified as:

  .. code-block:: xml

      <filter type="delayedgroup" bins="1 2 3 4 5 6" />

:energyfunction:
  ``energyfunction`` filters do not use the ``bins`` entry.  Instead
  they use ``energy`` and ``y``.


------------------
``<mesh>`` Element
------------------

If a structured mesh is desired as a filter for a tally, it must be specified in
a separate element with the tag name ``<mesh>``. This element has the following
attributes/sub-elements:

  :type:
    The type of structured mesh. The only valid option is "regular".

  :dimension:
    The number of mesh cells in each direction.

  :lower_left:
    The lower-left corner of the structured mesh. If only two coordinates are
    given, it is assumed that the mesh is an x-y mesh.

  :upper_right:
    The upper-right corner of the structured mesh. If only two coordinates are
    given, it is assumed that the mesh is an x-y mesh.

  :width:
    The width of mesh cells in each direction.

  .. note::
      One of ``<upper_right>`` or ``<width>`` must be specified, but not both
      (even if they are consistent with one another).

------------------------
``<derivative>`` Element
------------------------

OpenMC can take the first-order derivative of many tallies with respect to
material perturbations. It works by propagating a derivative through the
transport equation. Essentially, OpenMC keeps track of how each particle's
weight would change as materials are perturbed, and then accounts for that
weight change in the tallies. Note that this assumes material perturbations are
small enough not to change the distribution of fission sites. This element has
the following attributes/sub-elements:

  :id:
    A unique integer that can be used to identify the derivative.

  :variable:
    The independent variable of the derivative. Accepted options are "density",
    "nuclide_density", and "temperature". A "density" derivative will give the
    derivative with respect to the density of the material in [g / cm^3]. A
    "nuclide_density" derivative will give the derivative with respect to the
    density of a particular nuclide in units of [atom / b / cm].  A
    "temperature" derivative is with respect to a material temperature in units
    of [K].  The temperature derivative requires windowed multipole to be
    turned on.  Note also that the temperature derivative only accounts for
    resolved resonance Doppler broadening.  It does not account for thermal
    expansion, S(a, b) scattering, resonance scattering, or unresolved Doppler
    broadening.

  :material:
    The perturbed material. (Necessary for all derivative types)

  :nuclide:
    The perturbed nuclide. (Necessary only for "nuclide_density")

-----------------------------
``<assume_separate>`` Element
-----------------------------

In cases where the user needs to specify many different tallies each of which
are spatially separate, this tag can be used to cut down on some of the tally
overhead. The effect of assuming all tallies are spatially separate is that once
one tally is scored to, the same event is assumed not to score to any other
tallies. This element should be followed by "true" or "false".

  .. warning:: If used incorrectly, the assumption that all tallies are
               spatially separate can lead to incorrect results.

  *Default*: false
