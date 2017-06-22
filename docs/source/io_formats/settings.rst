.. _io_settings:

======================================
Settings Specification -- settings.xml
======================================

All simulation parameters and miscellaneous options are specified in the
settings.xml file.

---------------------
``<batches>`` Element
---------------------

The ``<batches>`` element indicates the total number of batches to execute,
where each batch corresponds to a tally realization. In a fixed source
calculation, each batch consists of a number of source particles. In an
eigenvalue calculation, each batch consists of one or many fission source
iterations (generations), where each generation itself consists of a number of
source neutrons.

  *Default*: None

----------------------------------
``<confidence_intervals>`` Element
----------------------------------

The ``<confidence_intervals>`` element has no attributes and has an accepted
value of "true" or "false". If set to "true", uncertainties on tally results
will be reported as the half-width of the 95% two-sided confidence interval. If
set to "false", uncertainties on tally results will be reported as the sample
standard deviation.

  *Default*: false

--------------------
``<cutoff>`` Element
--------------------

The ``<cutoff>`` element indicates two kinds of cutoffs. The first is the weight
cutoff used below which particles undergo Russian roulette. Surviving particles
are assigned a user-determined weight. Note that weight cutoffs and Russian
rouletting are not turned on by default. The second is the energy cutoff which
is used to kill particles under certain energy. The energy cutoff should not be
used unless you know particles under the energy are of no importance to results
you care. This element has the following attributes/sub-elements:

  :weight:
    The weight below which particles undergo Russian roulette.

    *Default*: 0.25

  :weight_avg:
    The weight that is assigned to particles that are not killed after Russian
    roulette.

    *Default*: 1.0

  :energy:
    The energy under which particles will be killed.

    *Default*: 0.0

-------------------------
``<energy_grid>`` Element
-------------------------

The ``<energy_grid>`` element determines the treatment of the energy grid during
a simulation. The valid options are "nuclide", "logarithm", and
"material-union". Setting this element to "nuclide" will cause OpenMC to use a
nuclide's energy grid when determining what points to interpolate between for
determining cross sections (i.e. non-unionized energy grid). Setting this
element to "logarithm" causes OpenMC to use a logarithmic mapping technique
described in LA-UR-14-24530_. Setting this element to "material-union" will
cause OpenMC to create energy grids that are unionized material-by-material and
use these grids when determining the energy-cross section pairs to interpolate
cross section values between.

  *Default*: logarithm

  .. note:: This element is not used in the multi-group :ref:`energy_mode`.

.. _LA-UR-14-24530: https://laws.lanl.gov/vhosts/mcnp.lanl.gov/pdf_files/la-ur-14-24530.pdf

.. _energy_mode:

-------------------------
``<energy_mode>`` Element
-------------------------

The ``<energy_mode>`` element tells OpenMC if the run-mode should be
continuous-energy or multi-group.  Options for entry are: ``continuous-energy``
or ``multi-group``.

  *Default*: continuous-energy

---------------------
``<entropy>`` Element
---------------------

The ``<entropy>`` element describes a mesh that is used for calculating Shannon
entropy. This mesh should cover all possible fissionable materials in the
problem. It has the following attributes/sub-elements:

  :dimension:
    The number of mesh cells in the x, y, and z directions, respectively.

    *Default*: If this tag is not present, the number of mesh cells is
     automatically determined by the code.

  :lower_left:
    The Cartesian coordinates of the lower-left corner of the mesh.

    *Default*: None

  :upper_right:
    The Cartesian coordinates of the upper-right corner of the mesh.

    *Default*: None

-----------------------------------
``<generations_per_batch>`` Element
-----------------------------------

The ``<generations_per_batch>`` element indicates the number of total fission
source iterations per batch for an eigenvalue calculation. This element is
ignored for all run modes other than "eigenvalue".

  *Default*: 1

----------------------
``<inactive>`` Element
----------------------

The ``<inactive>`` element indicates the number of inactive batches used in a
k-eigenvalue calculation. In general, the starting fission source iterations in
an eigenvalue calculation can not be used to contribute to tallies since the
fission source distribution and eigenvalue are generally not converged
immediately. This element is ignored for all run modes other than "eigenvalue".

  *Default*: 0

--------------------------
``<keff_trigger>`` Element
--------------------------

The ``<keff_trigger>`` element (ignored for all run modes other than
"eigenvalue".) specifies a precision trigger on the combined
:math:`k_{eff}`. The trigger is a convergence criterion on the uncertainty of
the estimated eigenvalue. It has the following attributes/sub-elements:

  :type:
    The type of precision trigger. Accepted options are "variance", "std_dev",
    and "rel_err".

    :variance:
      Variance of the batch mean :math:`\sigma^2`

    :std_dev:
      Standard deviation of the batch mean :math:`\sigma`

    :rel_err:
      Relative error of the batch mean :math:`\frac{\sigma}{\mu}`

    *Default*: None

  :threshold:
    The precision trigger's convergence criterion for the
    combined :math:`k_{eff}`.

    *Default*: None

.. note:: See section on the :ref:`trigger` for more information.


---------------------------
``<log_grid_bins>`` Element
---------------------------

The ``<log_grid_bins>`` element indicates the number of bins to use for the
logarithmic-mapped energy grid. Using more bins will result in energy grid
searches over a smaller range at the expense of more memory. The default is
based on the recommended value in LA-UR-14-24530_.

  *Default*: 8000

  .. note:: This element is not used in the multi-group :ref:`energy_mode`.

---------------------------
``<max_order>`` Element
---------------------------

The ``<max_order>`` element allows the user to set a maximum scattering order
to apply to every nuclide/material in the problem.  That is, if the data
library has :math:`P_3` data available, but ``<max_order>`` was set to ``1``,
then, OpenMC will only use up to the :math:`P_1` data.

  *Default*: Use the maximum order in the data library

  .. note:: This element is not used in the continuous-energy
    :ref:`energy_mode`.

-----------------------
``<no_reduce>`` Element
-----------------------

The ``<no_reduce>`` element has no attributes and has an accepted value of
"true" or "false". If set to "true", all user-defined tallies and global tallies
will not be reduced across processors in a parallel calculation. This means that
the accumulate score in one batch on a single processor is considered as an
independent realization for the tally random variable. For a problem with large
tally data, this option can significantly improve the parallel efficiency.

  *Default*: false

--------------------
``<output>`` Element
--------------------

The ``<output>`` element determines what output files should be written to disk
during the run. The sub-elements are described below, where "true" will write
out the file and "false" will not.

  :summary:
    Writes out an HDF5 summary file describing all of the user input files that
    were read in.

    *Default*: true

  :tallies:
    Write out an ASCII file of tally results.

    *Default*: true

  .. note:: The tally results will always be written to a binary/HDF5 state
            point file.

  :path:
    Absolute or relative path where all output files should be written to. The
    specified path must exist or else OpenMC will abort.

    *Default*: Current working directory

-----------------------
``<particles>`` Element
-----------------------

This element indicates the number of neutrons to simulate per fission source
iteration when a k-eigenvalue calculation is performed or the number of neutrons
per batch for a fixed source simulation.

  *Default*: None

---------------------
``<ptables>`` Element
---------------------

The ``<ptables>`` element determines whether probability tables should be used
in the unresolved resonance range if available. This element has no attributes
or sub-elements and can be set to either "false" or "true".

  *Default*: true

  .. note:: This element is not used in the multi-group :ref:`energy_mode`.

----------------------------------
``<resonance_scattering>`` Element
----------------------------------

The ``resonance_scattering`` element indicates to OpenMC that a method be used
to properly account for resonance elastic scattering (typically for nuclides
with Z > 40). This element can contain one or more of the following attributes
or sub-elements:

  :enable:
    Indicates whether a resonance elastic scattering method should be turned
    on. Accepts values of "true" or "false".

    *Default*: If the ``<resonance_scattering>`` element is present, "true".

  :method:

    Which resonance elastic scattering method is to be applied: "ares"
    (accelerated resonance elastic scattering), "dbrc" (Doppler broadening
    rejection correction), or "wcm" (weight correction method). Descriptions of
    each of these methods are documented here_.

    .. _here: http://dx.doi.org/10.1016/j.anucene.2014.01.017

    *Default*: "ares"

  :energy_min:
    The energy in eV above which the resonance elastic scattering method should
    be applied.

    *Default*: 0.01 eV

  :energy_max:
    The energy in eV below which the resonance elastic scattering method should
    be applied.

    *Default*: 1000.0 eV

  :nuclides:

    A list of nuclides to which the resonance elastic scattering method should
    be applied.

    *Default*: If ``<resonance_scattering>`` is present but the ``<nuclides>``
    sub-element is not given, the method is applied to all nuclides with 0 K
    elastic scattering data present.

  .. note:: If the ``resonance_scattering`` element is not given, the free gas,
            constant cross section scattering model, which has historically been
            used by Monte Carlo codes to sample target velocities, is used to
            treat the target motion of all nuclides.  If
            ``resonance_scattering`` is present, the constant cross section
            method is applied below ``energy_min`` and the target-at-rest
            (asymptotic) kernel is used above ``energy_max``.

  .. note:: This element is not used in the multi-group :ref:`energy_mode`.

----------------------
``<run_cmfd>`` Element
----------------------

The ``<run_cmfd>`` element indicates whether or not CMFD acceleration should be
turned on or off. This element has no attributes or sub-elements and can be set
to either "false" or "true".

  *Default*: false

----------------------
``<run_mode>`` Element
----------------------

The ``<run_mode>`` element indicates which run mode should be used when OpenMC
is executed. This element has no attributes or sub-elements and can be set to
"eigenvalue", "fixed source", "plot", "volume", or "particle restart".

  *Default*: None

------------------
``<seed>`` Element
------------------

The ``seed`` element is used to set the seed used for the linear congruential
pseudo-random number generator.

  *Default*: 1

--------------------
``<source>`` Element
--------------------

The ``source`` element gives information on an external source distribution to
be used either as the source for a fixed source calculation or the initial
source guess for criticality calculations. Multiple ``<source>`` elements may be
specified to define different source distributions. Each one takes the following
attributes/sub-elements:

  :strength:
    The strength of the source. If multiple sources are present, the source
    strength indicates the relative probability of choosing one source over the
    other.

    *Default*: 1.0

  :file:
    If this attribute is given, it indicates that the source is to be read from
    a binary source file whose path is given by the value of this element. Note,
    the number of source sites needs to be the same as the number of particles
    simulated in a fission source generation.

    *Default*: None

  :space:
    An element specifying the spatial distribution of source sites. This element
    has the following attributes:

    :type:
      The type of spatial distribution. Valid options are "box", "fission",
      "point", and "cartesian". A "box" spatial distribution has coordinates
      sampled uniformly in a parallelepiped. A "fission" spatial distribution
      samples locations from a "box" distribution but only locations in
      fissionable materials are accepted. A "point" spatial distribution has
      coordinates specified by a triplet. An "cartesian" spatial distribution
      specifies independent distributions of x-, y-, and z-coordinates.

      *Default*: None

    :parameters:
      For a "box" or "fission" spatial distribution, ``parameters`` should be
      given as six real numbers, the first three of which specify the lower-left
      corner of a parallelepiped and the last three of which specify the
      upper-right corner. Source sites are sampled uniformly through that
      parallelepiped.

      For a "point" spatial distribution, ``parameters`` should be given as
      three real numbers which specify the (x,y,z) location of an isotropic
      point source.

      For an "cartesian" distribution, no parameters are specified. Instead,
      the ``x``, ``y``, and ``z`` elements must be specified.

      *Default*: None

    :x:
      For an "cartesian" distribution, this element specifies the distribution
      of x-coordinates. The necessary sub-elements/attributes are those of a
      univariate probability distribution (see the description in
      :ref:`univariate`).

    :y:
      For an "cartesian" distribution, this element specifies the distribution
      of y-coordinates. The necessary sub-elements/attributes are those of a
      univariate probability distribution (see the description in
      :ref:`univariate`).

    :z:
      For an "cartesian" distribution, this element specifies the distribution
      of z-coordinates. The necessary sub-elements/attributes are those of a
      univariate probability distribution (see the description in
      :ref:`univariate`).

  :angle:
    An element specifying the angular distribution of source sites. This element
    has the following attributes:

    :type:
      The type of angular distribution. Valid options are "isotropic",
      "monodirectional", and "mu-phi". The angle of the particle emitted from a
      source site is isotropic if the "isotropic" option is given. The angle of
      the particle emitted from a source site is the direction specified in the
      ``reference_uvw`` element/attribute if "monodirectional" option is
      given. The "mu-phi" option produces directions with the cosine of the
      polar angle and the azimuthal angle explicitly specified.

      *Default*: isotropic

    :reference_uvw:
      The direction from which the polar angle is measured. Represented by the
      x-, y-, and z-components of a unit vector. For a monodirectional
      distribution, this defines the direction of all sampled particles.

    :mu:
      An element specifying the distribution of the cosine of the polar
      angle. Only relevant when the type is "mu-phi". The necessary
      sub-elements/attributes are those of a univariate probability distribution
      (see the description in :ref:`univariate`).

    :phi:
      An element specifying the distribution of the azimuthal angle. Only
      relevant when the type is "mu-phi". The necessary sub-elements/attributes
      are those of a univariate probability distribution (see the description in
      :ref:`univariate`).

  :energy:
    An element specifying the energy distribution of source sites. The necessary
    sub-elements/attributes are those of a univariate probability distribution
    (see the description in :ref:`univariate`).

    *Default*: Watt spectrum with :math:`a` = 0.988 MeV and :math:`b` =
    2.249 MeV :sup:`-1`

  :write_initial:
    An element specifying whether to write out the initial source bank used at
    the beginning of the first batch. The output file is named
    "initial_source.h5"

    *Default*: false

.. _univariate:

Univariate Probability Distributions
++++++++++++++++++++++++++++++++++++

Various components of a source distribution involve probability distributions of
a single random variable, e.g. the distribution of the energy, the distribution
of the polar angle, and the distribution of x-coordinates. Each of these
components supports the same syntax with an element whose tag signifies the
variable and whose sub-elements/attributes are as follows:

:type:
  The type of the distribution. Valid options are "uniform", "discrete",
  "tabular", "maxwell", and "watt". The "uniform" option produces variates
  sampled from a uniform distribution over a finite interval. The "discrete"
  option produces random variates that can assume a finite number of values
  (i.e., a distribution characterized by a probability mass function). The
  "tabular" option produces random variates sampled from a tabulated
  distribution where the density function is either a histogram or
  linearly-interpolated between tabulated points. The "watt" option produces
  random variates is sampled from a Watt fission spectrum (only used for
  energies). The "maxwell" option produce variates sampled from a Maxwell
  fission spectrum (only used for energies).

  *Default*: None

:parameters:
  For a "uniform" distribution, ``parameters`` should be given as two real
  numbers :math:`a` and :math:`b` that define the interval :math:`[a,b]` over
  which random variates are sampled.

  For a "discrete" or "tabular" distribution, ``parameters`` provides the
  :math:`(x,p)` pairs defining the discrete/tabular distribution. All :math:`x`
  points are given first followed by corresponding :math:`p` points.

  For a "watt" distribution, ``parameters`` should be given as two real numbers
  :math:`a` and :math:`b` that parameterize the distribution :math:`p(x) dx = c
  e^{-x/a} \sinh \sqrt{b \, x} dx`.

  For a "maxwell" distribution, ``parameters`` should be given as one real
  number :math:`a` that parameterizes the distribution :math:`p(x) dx = c x
  e^{-x/a} dx`.

  .. note:: The above format should be used even when using the multi-group
            :ref:`energy_mode`.
:interpolation:
  For a "tabular" distribution, ``interpolation`` can be set to "histogram" or
  "linear-linear" thereby specifying how tabular points are to be interpolated.

  *Default*: histogram

-------------------------
``<state_point>`` Element
-------------------------

The ``<state_point>`` element indicates at what batches a state point file
should be written. A state point file can be used to restart a run or to get
tally results at any batch. The default behavior when using this tag is to
write out the source bank in the state_point file. This behavior can be
customized by using the ``<source_point>`` element. This element has the
following attributes/sub-elements:

  :batches:
    A list of integers separated by spaces indicating at what batches a state
    point file should be written.

    *Default*: Last batch only

--------------------------
``<source_point>`` Element
--------------------------

The ``<source_point>`` element indicates at what batches the source bank
should be written. The source bank can be either written out within a state
point file or separately in a source point file. This element has the following
attributes/sub-elements:

  :batches:
    A list of integers separated by spaces indicating at what batches a state
    point file should be written. It should be noted that if the ``separate``
    attribute is not set to "true", this list must be a subset of state point
    batches.

    *Default*: Last batch only

  :separate:
    If this element is set to "true", a separate binary source point file will
    be written. Otherwise, the source sites will be written in the state point
    directly.

    *Default*: false

  :write:
    If this element is set to "false", source sites are not written
    to the state point or source point file. This can substantially reduce the
    size of state points if large numbers of particles per batch are used.

    *Default*: true

  :overwrite_latest:
    If this element is set to "true", a source point file containing
    the source bank will be written out to a separate file named
    ``source.binary`` or ``source.h5`` depending on if HDF5 is enabled.
    This file will be overwritten at every single batch so that the latest
    source bank will be available. It should be noted that a user can set both
    this element to "true" and specify batches to write a permanent source bank.

    *Default*: false

------------------------------
``<survival_biasing>`` Element
------------------------------

The ``<survival_biasing>`` element has no attributes and has an accepted value
of "true" or "false". If set to "true", this option will enable the use of
survival biasing, otherwise known as implicit capture or absorption.

  *Default*: false

.. _tabular_legendre:

---------------------------------
``<tabular_legendre>`` Element
---------------------------------

The optional ``<tabular_legendre>`` element specifies how the multi-group
Legendre scattering kernel is represented if encountered in a multi-group
problem.  Specifically, the options are to either convert the Legendre
expansion to a tabular representation or leave it as a set of Legendre
coefficients. Converting to a tabular representation will cost memory but can
allow for a decrease in runtime compared to leaving as a set of Legendre
coefficients. This element has the following attributes/sub-elements:

  :enable:
    This attribute/sub-element denotes whether or not the conversion of a
    Legendre scattering expansion to the tabular format should be performed or
    not. A value of “true” means the conversion should be performed, “false”
    means it will not.

    *Default*: true

  :num_points:
    If the conversion is to take place the number of tabular points is
    required. This attribute/sub-element allows the user to set the desired
    number of points.

    *Default*: 33

  .. note:: This element is only used in the multi-group :ref:`energy_mode`.

.. _temperature_default:

---------------------------------
``<temperature_default>`` Element
---------------------------------

The ``<temperature_default>`` element specifies a default temperature in Kelvin
that is to be applied to cells in the absence of an explicit cell temperature or
a material default temperature.

  *Default*: 293.6 K

.. _temperature_method:

--------------------------------
``<temperature_method>`` Element
--------------------------------

The ``<temperature_method>`` element has an accepted value of "nearest" or
"interpolation". A value of "nearest" indicates that for each
cell, the nearest temperature at which cross sections are given is to be
applied, within a given tolerance (see :ref:`temperature_tolerance`). A value of
"interpolation" indicates that cross sections are to be linear-linear
interpolated between temperatures at which nuclear data are present (see
:ref:`temperature_treatment`).

  *Default*: "nearest"

.. _temperature_multipole:

-----------------------------------
``<temperature_multipole>`` Element
-----------------------------------

The ``<temperature_multipole>`` element toggles the windowed multipole
capability on or off. If this element is set to "True" and the relevant data is
available, OpenMC will use the windowed multipole method to evaluate and Doppler
broaden cross sections in the resolved resonance range.  This override other
methods like "nearest" and "interpolation" in the resolved resonance range.

  *Default*: False

-------------------------------
``<temperature_range>`` Element
-------------------------------

The ``<temperature_range>`` element specifies a minimum and maximum temperature
in Kelvin above and below which cross sections should be loaded for all nuclides
and thermal scattering tables. This can be used for multi-physics simulations
where the temperatures might change from one iteration to the next.

  *Default*: None

.. _temperature_tolerance:

-----------------------------------
``<temperature_tolerance>`` Element
-----------------------------------

The ``<temperature_tolerance>`` element specifies a tolerance in Kelvin that is
to be applied when the "nearest" temperature method is used. For example, if a
cell temperature is 340 K and the tolerance is 15 K, then the closest
temperature in the range of 325 K to 355 K will be used to evaluate cross
sections.

  *Default*: 10 K

---------------------
``<threads>`` Element
---------------------

The ``<threads>`` element indicates the number of OpenMP threads to be used for
a simulation. It has no attributes and accepts a positive integer value.

  *Default*: None (Determined by environment variable :envvar:`OMP_NUM_THREADS`)

.. _trace:

-------------------
``<trace>`` Element
-------------------

The ``<trace>`` element can be used to print out detailed information about a
single particle during a simulation. This element should be followed by three
integers: the batch number, generation number, and particle number.

  *Default*: None

.. _track:

-------------------
``<track>`` Element
-------------------

The ``<track>`` element specifies particles for which OpenMC will output binary
files describing particle position at every step of its transport. This element
should be followed by triplets of integers.  Each triplet describes one
particle. The integers in each triplet specify the batch number, generation
number, and particle number, respectively.

  *Default*: None

.. _trigger:

-------------------------
``<trigger>`` Element
-------------------------

OpenMC includes tally precision triggers which allow the user to define
uncertainty thresholds on :math:`k_{eff}` in the ``<keff_trigger>`` subelement
of ``settings.xml``, and/or tallies in ``tallies.xml``. When using triggers,
OpenMC will run until it completes as many batches as defined by ``<batches>``.
At this point, the uncertainties on all tallied values are computed and compared
with their corresponding trigger thresholds. If any triggers have not been met,
OpenMC will continue until either all trigger thresholds have been satisfied or
``<max_batches>`` has been reached.

The ``<trigger>`` element provides an active "toggle switch" for tally
precision trigger(s), the maximum number of batches and the batch interval. It
has the following attributes/sub-elements:

  :active:
    This determines whether or not to use trigger(s). Trigger(s) are used when
    this tag is set to "true".

  :max_batches:
    This describes the maximum number of batches allowed when using trigger(s).

    .. note:: When max_batches is set, the number of ``batches`` shown in the
              ``<batches>`` element represents minimum number of batches to
              simulate when using the trigger(s).

  :batch_interval:
    This tag describes the number of  batches in between convergence checks.
    OpenMC will check if the trigger has been reached at each batch defined
    by ``batch_interval`` after the minimum number of batches is reached.

    .. note:: If this tag is not present, the ``batch_interval`` is predicted
              dynamically by OpenMC for each convergence check. The predictive
              model assumes no correlation between fission sources
              distributions from batch-to-batch. This assumption is reasonable
              for fixed source and small criticality calculations, but is very
              optimistic for highly coupled full-core reactor problems.


------------------------
``<uniform_fs>`` Element
------------------------

The ``<uniform_fs>`` element describes a mesh that is used for re-weighting
source sites at every generation based on the uniform fission site methodology
described in Kelly et al., "MC21 Analysis of the Nuclear Energy Agency Monte
Carlo Performance Benchmark Problem," Proceedings of *Physor 2012*, Knoxville,
TN (2012). This mesh should cover all possible fissionable materials in the
problem. It has the following attributes/sub-elements:

  :dimension:
    The number of mesh cells in the x, y, and z directions, respectively.

    *Default*: None

  :lower_left:
    The Cartesian coordinates of the lower-left corner of the mesh.

    *Default*: None

  :upper_right:
    The Cartesian coordinates of the upper-right corner of the mesh.

    *Default*: None

.. _verbosity:

-----------------------
``<verbosity>`` Element
-----------------------

The ``<verbosity>`` element tells the code how much information to display to
the standard output. A higher verbosity corresponds to more information being
displayed. The text of this element should be an integer between between 1
and 10. The verbosity levels are defined as follows:

  :1: don't display any output
  :2: only show OpenMC logo
  :3: all of the above + headers
  :4: all of the above + results
  :5: all of the above + file I/O
  :6: all of the above + timing statistics and initialization messages
  :7: all of the above + :math:`k` by generation
  :9: all of the above + indicate when each particle starts
  :10: all of the above + event information

  *Default*: 7

-------------------------------------
``<create_fission_neutrons>`` Element
-------------------------------------

The ``<create_fission_neutrons>`` element indicates whether fission neutrons
should be created or not.  If this element is set to "true", fission neutrons
will be created; otherwise the fission is treated as capture and no fission
neutron will be created. Note that this option is only applied to fixed source
calculation. For eigenvalue calculation, fission will always be treated as real
fission.

  *Default*: true


-------------------------
``<volume_calc>`` Element
-------------------------

The ``<volume_calc>`` element indicates that a stochastic volume calculation
should be run at the beginning of the simulation. This element has the following
sub-elements/attributes:

  :cells:
    The unique IDs of cells for which the volume should be estimated.

    *Default*: None

  :samples:
    The number of samples used to estimate volumes.

    *Default*: None

  :lower_left:
     The lower-left Cartesian coordinates of a bounding box that is used to
     sample points within.

     *Default*: None

  :upper_right:
     The upper-right Cartesian coordinates of a bounding box that is used to
     sample points within.

     *Default*: None
