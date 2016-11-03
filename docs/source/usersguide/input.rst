.. _usersguide_input:

=======================
Writing XML Input Files
=======================

Unlike many other Monte Carlo codes which use an arbitrary-format ASCII file
with "cards" to specify a particular geometry, materials, and associated run
settings, the input files for OpenMC are structured in a set of XML_ files. XML,
which stands for eXtensible Markup Language, is a simple format that allows data
to be exchanged efficiently between different programs and interfaces.

Anyone who has ever seen webpages written in HTML will be familiar with the
structure of XML whereby "tags" enclosed in angle brackets denote that a
particular piece of data will follow. Let us examine the follow example:

.. code-block:: xml

    <person>
      <firstname>John</firstname>
      <lastname>Smith</lastname>
      <age>27</age>
      <occupation>Health Physicist</occupation>
    </person>

Here we see that the first tag indicates that the following data will describe a
person. The nested tags *firstname*, *lastname*, *age*, and *occupation*
indicate characteristics about the person being described.

In much the same way, OpenMC input uses XML tags to describe the geometry, the
materials, and settings for a Monte Carlo simulation.

.. _XML: http://www.w3.org/XML/

-----------------
Overview of Files
-----------------

To assemble a complete model for OpenMC, one needs to create separate XML files
for the geometry, materials, and settings. Additionally, there are three optional
input files. The first is a tallies XML file that specifies physical quantities
to be tallied. The second is a plots XML file that specifies regions of geometry
which should be plotted. The third is a CMFD XML file that specifies coarse mesh
acceleration geometry and execution parameters. OpenMC expects that these
files are called:

* ``geometry.xml``
* ``materials.xml``
* ``settings.xml``
* ``tallies.xml``
* ``plots.xml``
* ``cmfd.xml``

--------------------
Validating XML Files
--------------------

Input files can be checked before executing OpenMC using the
``openmc-validate-xml`` script which is installed alongside the Python API. Two
command line arguments can be set when running ``openmc-validate-xml``:

* ``-i``, ``--input-path`` - Location of OpenMC input files.
  *Default*: current working directory
* ``-r``, ``--relaxng-path`` - Location of OpenMC RelaxNG files.
  *Default*: None

If the RelaxNG path is not set, the script will search for these files because
it expects that the user is either running the script located in the install
directory ``bin`` folder or in ``src/utils``. Once executed, it will match
OpenMC XML files with their RelaxNG schema and check if they are valid.  Below
is a table of the messages that will be printed after each file is checked.

========================  ===================================
Message                   Description
========================  ===================================
[XML ERROR]               Cannot parse XML file.
[NO RELAXNG FOUND]        No RelaxNG file found for XML file.
[NOT VALID]               XML file does not match RelaxNG.
[VALID]                   XML file matches RelaxNG.
========================  ===================================

As an example, if OpenMC is installed in the directory ``/opt/openmc/`` and the
current working directory is where OpenMC XML input files are located, they can
be validated using the following command:

.. code-block:: bash

   /opt/openmc/bin/openmc-validate-xml

--------------------------------------
Settings Specification -- settings.xml
--------------------------------------

All simulation parameters and miscellaneous options are specified in the
settings.xml file.

``<confidence_intervals>`` Element
----------------------------------

The ``<confidence_intervals>`` element has no attributes and has an accepted
value of "true" or "false". If set to "true", uncertainties on tally results
will be reported as the half-width of the 95% two-sided confidence interval. If
set to "false", uncertainties on tally results will be reported as the sample
standard deviation.

  *Default*: false

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

.. _eigenvalue:

``<eigenvalue>`` Element
------------------------

The ``<eigenvalue>`` element indicates that a :math:`k`-eigenvalue calculation
should be performed. It has the following attributes/sub-elements:

  :batches:
    The total number of batches, where each batch corresponds to multiple
    fission source iterations. Batching is done to eliminate correlation between
    realizations of random variables.

    *Default*: None

  :generations_per_batch:
    The number of total fission source iterations per batch.

    *Default*: 1

  :inactive:
    The number of inactive batches. In general, the starting cycles in a
    criticality calculation can not be used to contribute to tallies since the
    fission source distribution and eigenvalue are generally not converged
    immediately.

    *Default*: None

  :particles:
    The number of neutrons to simulate per fission source iteration.

    *Default*: None

  :keff_trigger:
    This tag specifies a precision trigger on the combined :math:`k_{eff}`. The
    trigger is a convergence criterion on the uncertainty of the estimated
    eigenvalue. It has the following attributes/sub-elements:

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

``<energy_mode>`` Element
-------------------------

The ``<energy_mode>`` element tells OpenMC if the run-mode should be
continuous-energy or multi-group.  Options for entry are: ``continuous-energy``
or ``multi-group``.

  *Default*: continuous-energy

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

``<fixed_source>`` Element
--------------------------

The ``<fixed_source>`` element indicates that a fixed source calculation should
be performed. It has the following attributes/sub-elements:

  :batches:
    The total number of batches. For fixed source calculations, each batch
    represents a realization of random variables for tallies.

    *Default*: None

  :particles:
    The number of particles to simulate per batch.

    *Default*: None

``<log_grid_bins>`` Element
---------------------------

The ``<log_grid_bins>`` element indicates the number of bins to use for the
logarithmic-mapped energy grid. Using more bins will result in energy grid
searches over a smaller range at the expense of more memory. The default is
based on the recommended value in LA-UR-14-24530_.

  *Default*: 8000

  .. note:: This element is not used in the multi-group :ref:`energy_mode`.

``<max_order>`` Element
---------------------------

The ``<max_order>`` element allows the user to set a maximum scattering order
to apply to every nuclide/material in the problem.  That is, if the data
library has :math:`P_3` data available, but ``<max_order>`` was set to ``1``,
then, OpenMC will only use up to the :math:`P_1` data.

  *Default*: Use the maximum order in the data library

  .. note:: This element is not used in the continuous-energy
    :ref:`energy_mode`.

``<no_reduce>`` Element
-----------------------

The ``<no_reduce>`` element has no attributes and has an accepted value of
"true" or "false". If set to "true", all user-defined tallies and global tallies
will not be reduced across processors in a parallel calculation. This means that
the accumulate score in one batch on a single processor is considered as an
independent realization for the tally random variable. For a problem with large
tally data, this option can significantly improve the parallel efficiency.

  *Default*: false

``<output>`` Element
--------------------

The ``<output>`` element determines what output files should be written to disk
during the run. The sub-elements are described below, where "true" will write
out the file and "false" will not.

  :cross_sections:
    Writes out an ASCII summary file of the cross sections that were read in.

    *Default*: false

  :summary:
    Writes out an HDF5 summary file describing all of the user input files that
    were read in.

    *Default*: true

  :tallies:
    Write out an ASCII file of tally results.

    *Default*: true

  .. note:: The tally results will always be written to a binary/HDF5 state
            point file.

``<output_path>`` Element
-------------------------

The ``<output_path>`` element specifies an absolute or relative path where all
output files should be written to. The specified path must exist or else OpenMC
will abort.

  *Default*: Current working directory

``<ptables>`` Element
---------------------

The ``<ptables>`` element determines whether probability tables should be used
in the unresolved resonance range if available. This element has no attributes
or sub-elements and can be set to either "false" or "true".

  *Default*: true

  .. note:: This element is not used in the multi-group :ref:`energy_mode`.

``<resonance_scattering>`` Element
----------------------------------

The ``resonance_scattering`` element can contain one or more of the following
attributes or sub-elements:

  :scatterer:
    An element with attributes/sub-elements called ``nuclide``, ``method``,
    ``E_min``, and ``E_max``. The ``nuclide`` attribute is the name, as given
    by the ``name`` attribute within the ``nuclide`` sub-element of the
    ``material`` element in ``materials.xml``, of the nuclide to which a
    resonance scattering treatment is to be applied.
    The ``method`` attribute gives the type of resonance scattering treatment
    that is to be applied to the ``nuclide``.  Acceptable inputs - none of
    which are case-sensitive - for the ``method`` attribute are ``ARES``,
    ``CXS``, ``WCM``, and ``DBRC``.  Descriptions of each of these methods
    are documented here_.  The ``E_min`` attribute gives the minimum energy
    above which the ``method`` is applied.  The ``E_max`` attribute gives the
    maximum energy below which the ``method`` is applied.  One example would
    be as follows:

    .. _here: http://dx.doi.org/10.1016/j.anucene.2014.01.017

    .. code-block:: xml

        <resonance_scattering>
          <scatterer>
            <nuclide>U-238</nuclide>
            <method>ARES</method>
            <E_min>5.0e-6</E_min>
            <E_max>40.0e-6</E_max>
         </scatterer>
         <scatterer>
            <nuclide>Pu-239</nuclide>
            <method>dbrc</method>
            <E_min>0.01e-6</E_min>
            <E_max>210.0e-6</E_max>
          </scatterer>
        </resonance_scattering>

    .. note:: If the ``resonance_scattering`` element is not given, the free gas,
              constant cross section (``cxs``) scattering model, which has
              historically been used by Monte Carlo codes to sample target
              velocities, is used to treat the target motion of all nuclides.  If
              ``resonance_scattering`` is present, the ``cxs`` method is applied
              below ``E_min`` and the target-at-rest (asymptotic) kernel is used
              above ``E_max``.  An arbitrary number of ``scatterer`` elements may
              be specified, each corresponding to a single nuclide at a single
              temperature.

    *Defaults*: None (scatterer), ARES (method), 0.01 eV (E_min), 1.0 keV (E_max)

  .. note:: This element is not used in the multi-group :ref:`energy_mode`.

``<run_cmfd>`` Element
----------------------

The ``<run_cmfd>`` element indicates whether or not CMFD acceleration should be
turned on or off. This element has no attributes or sub-elements and can be set
to either "false" or "true".

  *Defualt*: false

``<seed>`` Element
------------------

The ``seed`` element is used to set the seed used for the linear congruential
pseudo-random number generator.

  *Default*: 1

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

  :interval:
    A single integer :math:`n` indicating that a state point should be written
    every :math:`n` batches. This option can be given in lieu of listing
    batches explicitly.

    *Default*: None

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

  :interval:
    A single integer :math:`n` indicating that a state point should be written
    every :math:`n` batches. This option can be given in lieu of listing batches
    explicitly. It should be noted that if the ``separate`` attribute is not set
    to "true", this value should produce a list of batches that is a subset of
    state point batches.

    *Default*: None

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

``<survival_biasing>`` Element
------------------------------

The ``<survival_biasing>`` element has no attributes and has an accepted value
of "true" or "false". If set to "true", this option will enable the use of
survival biasing, otherwise known as implicit capture or absorption.

  *Default*: false

.. _tabular_legendre:

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

``<temperature_default>`` Element
---------------------------------

The ``<temperature_default>`` element specifies a default temperature in Kelvin
that is to be applied to cells in the absence of an explicit cell temperature or
a material default temperature.

  *Default*: 293.6 K

.. _temperature_method:

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

``<temperature_multipole>`` Element
-----------------------------------

The ``<temperature_multipole>`` element toggles the windowed multipole
capability on or off. If this element is set to "True" and the relevant data is
available, OpenMC will use the windowed multipole method to evaluate and Doppler
broaden cross sections in the resolved resonance range.  This override other
methods like "nearest" and "interpolation" in the resolved resonance range.

  *Default*: False

.. _temperature_tolerance:

``<temperature_tolerance>`` Element
-----------------------------------

The ``<temperature_tolerance>`` element specifies a tolerance in Kelvin that is
to be applied when the "nearest" temperature method is used. For example, if a
cell temperature is 340 K and the tolerance is 15 K, then the closest
temperature in the range of 325 K to 355 K will be used to evaluate cross
sections.

  *Default*: 10 K

``<threads>`` Element
---------------------

The ``<threads>`` element indicates the number of OpenMP threads to be used for
a simulation. It has no attributes and accepts a positive integer value.

  *Default*: None (Determined by environment variable :envvar:`OMP_NUM_THREADS`)

.. _trace:

``<trace>`` Element
-------------------

The ``<trace>`` element can be used to print out detailed information about a
single particle during a simulation. This element should be followed by three
integers: the batch number, generation number, and particle number.

  *Default*: None

.. _track:

``<track>`` Element
-------------------

The ``<track>`` element specifies particles for which OpenMC will output binary
files describing particle position at every step of its transport. This element
should be followed by triplets of integers.  Each triplet describes one
particle. The integers in each triplet specify the batch number, generation
number, and particle number, respectively.

  *Default*: None

.. _trigger:

``<trigger>`` Element
-------------------------

OpenMC includes tally precision triggers which allow the user to define
uncertainty thresholds on :math:`k_{eff}` in the ``<eigenvalue>`` subelement of
``settings.xml``, and/or tallies in ``tallies.xml``. When using triggers,
OpenMC will run until it completes as many batches as defined by ``<batches>``.
At this point, the uncertainties on all tallied values are computed and
compared with their corresponding trigger thresholds. If any triggers have not
been met, OpenMC will continue until either all trigger thresholds have been
satisfied or ``<max_batches>`` has been reached.

The ``<trigger>`` element provides an active "toggle switch" for tally
precision trigger(s), the maximum number of batches and the batch interval. It
has the following attributes/sub-elements:

  :active:
    This determines whether or not to use trigger(s). Trigger(s) are used when
    this tag is set to "true".

  :max_batches:
    This describes the maximum number of batches allowed when using trigger(s).

    .. note:: When max_batches is set, the number of ``batches`` shown in
              ``<eigenvalue>`` element represents minimum number of batches to
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

``<verbosity>`` Element
-----------------------

The ``<verbosity>`` element tells the code how much information to display to
the standard output. A higher verbosity corresponds to more information being
displayed. This element takes the following attributes:

  :value:
    The specified verbosity between 1 and 10.

    *Default*: 5

``<create_fission_neutrons>`` Element
-------------------------------------

The ``<create_fission_neutrons>`` element indicates whether fission neutrons
should be created or not.  If this element is set to "true", fission neutrons
will be created; otherwise the fission is treated as capture and no fission
neutron will be created. Note that this option is only applied to fixed source
calculation. For eigenvalue calculation, fission will always be treated as real
fission.

  *Default*: true


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

--------------------------------------
Geometry Specification -- geometry.xml
--------------------------------------

The geometry in OpenMC is described using `constructive solid geometry`_ (CSG),
also sometimes referred to as combinatorial geometry. CSG allows a user to
create complex objects using Boolean operators on a set of simpler surfaces. In
the geometry model, each unique volume is defined by its bounding surfaces. In
OpenMC, most `quadratic surfaces`_ can be modeled and used as bounding surfaces.

Every geometry.xml must have an XML declaration at the beginning of the file and
a root element named geometry. Within the root element the user can define any
number of cells, surfaces, and lattices. Let us look at the following example:

.. code-block:: xml

    <?xml version="1.0"?>
    <geometry>
      <!-- This is a comment -->

      <surface>
        <id>1</id>
        <type>sphere</type>
        <coeffs>0.0 0.0 0.0 5.0</coeffs>
        <boundary>vacuum</boundary>
      <surface>

      <cell>
        <id>1</id>
        <universe>0</universe>
        <material>1</material>
        <region>-1</region>
      </cell>
    </geometry>

At the beginning of this file is a comment, denoted by a tag starting with
``<!--`` and ending with ``-->``. Comments, as well as any other type of input,
may span multiple lines. One convenient feature of the XML input format is that
sub-elements of the ``cell`` and ``surface`` elements can also be equivalently
expressed of attributes of the original element, e.g. the geometry file above
could be written as:

.. code-block:: xml

    <?xml version="1.0"?>
    <geometry>
      <!-- This is a comment -->

      <surface id="1" type="sphere" coeffs="0.0 0.0 0.0 5.0" boundary="vacuum" />
      <cell id="1" universe="0" material="1" region="-1" />

    </geometry>

.. _surface_element:

``<surface>`` Element
---------------------

Each ``<surface>`` element can have the following attributes or sub-elements:

  :id:
    A unique integer that can be used to identify the surface.

    *Default*: None

  :name:
    An optional string name to identify the surface in summary output
    files. This string is limited to 52 characters for formatting purposes.

    *Default*: ""

  :type:
    The type of the surfaces. This can be "x-plane", "y-plane", "z-plane",
    "plane", "x-cylinder", "y-cylinder", "z-cylinder", "sphere", "x-cone",
    "y-cone", "z-cone", or "quadric".

    *Default*: None

  :coeffs:
    The corresponding coefficients for the given type of surface. See below for
    a list a what coefficients to specify for a given surface

    *Default*: None

  :boundary:
     The boundary condition for the surface. This can be "transmission",
     "vacuum", "reflective", or "periodic". Periodic boundary conditions can
     only be applied to x-, y-, and z-planes. Only axis-aligned periodicity is
     supported, i.e., x-planes can only be paired with x-planes. Specify which
     planes are periodic and the code will automatically identify which planes
     are paired together.

    *Default*: "transmission"

  :periodic_surface_id:
     If a periodic boundary condition is applied, this attribute identifies the
     ``id`` of the corresponding periodic sufrace.

The following quadratic surfaces can be modeled:

  :x-plane:
    A plane perpendicular to the x axis, i.e. a surface of the form :math:`x -
    x_0 = 0`. The coefficients specified are ":math:`x_0`".

  :y-plane:
    A plane perpendicular to the y axis, i.e. a surface of the form :math:`y -
    y_0 = 0`. The coefficients specified are ":math:`y_0`".

  :z-plane:
    A plane perpendicular to the z axis, i.e. a surface of the form :math:`z -
    z_0 = 0`. The coefficients specified are ":math:`z_0`".

  :plane:
    An arbitrary plane of the form :math:`Ax + By + Cz = D`. The coefficients
    specified are ":math:`A \: B \: C \: D`".

  :x-cylinder:
    An infinite cylinder whose length is parallel to the x-axis. This is a
    quadratic surface of the form :math:`(y - y_0)^2 + (z - z_0)^2 = R^2`. The
    coefficients specified are ":math:`y_0 \: z_0 \: R`".

  :y-cylinder:
    An infinite cylinder whose length is parallel to the y-axis. This is a
    quadratic surface of the form :math:`(x - x_0)^2 + (z - z_0)^2 = R^2`. The
    coefficients specified are ":math:`x_0 \: z_0 \: R`".

  :z-cylinder:
    An infinite cylinder whose length is parallel to the z-axis. This is a
    quadratic surface of the form :math:`(x - x_0)^2 + (y - y_0)^2 = R^2`. The
    coefficients specified are ":math:`x_0 \: y_0 \: R`".

  :sphere:
    A sphere of the form :math:`(x - x_0)^2 + (y - y_0)^2 + (z - z_0)^2 =
    R^2`. The coefficients specified are ":math:`x_0 \: y_0 \: z_0 \: R`".

  :x-cone:
    A cone parallel to the x-axis of the form :math:`(y - y_0)^2 + (z - z_0)^2 =
    R^2 (x - x_0)^2`. The coefficients specified are ":math:`x_0 \: y_0 \: z_0
    \: R^2`".

  :y-cone:
    A cone parallel to the y-axis of the form :math:`(x - x_0)^2 + (z - z_0)^2 =
    R^2 (y - y_0)^2`. The coefficients specified are ":math:`x_0 \: y_0 \: z_0
    \: R^2`".

  :z-cone:
    A cone parallel to the x-axis of the form :math:`(x - x_0)^2 + (y - y_0)^2 =
    R^2 (z - z_0)^2`. The coefficients specified are ":math:`x_0 \: y_0 \: z_0
    \: R^2`".

  :quadric:
     A general quadric surface of the form :math:`Ax^2 + By^2 + Cz^2 + Dxy +
     Eyz + Fxz + Gx + Hy + Jz + K = 0` The coefficients specified are ":math:`A
     \: B \: C \: D \: E \: F \: G \: H \: J \: K`".


``<cell>`` Element
------------------

Each ``<cell>`` element can have the following attributes or sub-elements:

  :id:
    A unique integer that can be used to identify the cell.

    *Default*: None

  :name:
    An optional string name to identify the cell in summary output files.
    This string is limmited to 52 characters for formatting purposes.

    *Default*: ""

  :universe:
    The ``id`` of the universe that this cell is contained in.

    *Default*: 0

  :fill:
    The ``id`` of the universe that fills this cell.

    .. note:: If a fill is specified, no material should be given.

    *Default*: None

  :material:
    The ``id`` of the material that this cell contains. If the cell should
    contain no material, this can also be set to "void". A list of materials
    can be specified for the "distributed material" feature. This will give each
    unique instance of the cell its own material.

    .. note:: If a material is specified, no fill should be given.

    *Default*: None

  :region:
    A Boolean expression of half-spaces that defines the spatial region which
    the cell occupies. Each half-space is identified by the unique ID of the
    surface prefixed by `-` or `+` to indicate that it is the negative or
    positive half-space, respectively. The `+` sign for a positive half-space
    can be omitted. Valid Boolean operators are parentheses, union `|`,
    complement `~`, and intersection. Intersection is implicit and indicated by
    the presence of whitespace. The order of operator precedence is parentheses,
    complement, intersection, and then union.

    As an example, the following code gives a cell that is the union of the
    negative half-space of surface 3 and the complement of the intersection of
    the positive half-space of surface 5 and the negative half-space of surface
    2:

    .. code-block:: xml

        <cell id="1" material="1" region="-3 | ~(5 -2)" />

    .. note:: The ``region`` attribute/element can be omitted to make a cell
              fill its entire universe.

    *Default*: A region filling all space.

  :temperature:
    The temperature of the cell in Kelvin. If windowed-multipole data is
    avalable, this temperature will be used to Doppler broaden some cross
    sections in the resolved resonance region. A list of temperatures can be
    specified for the "distributed temperature" feature. This will give each
    unique instance of the cell its own temperature.

    *Default*: If a material default temperature is supplied, it is used. In the
    absence of a material default temperature, the :ref:`global default
    temperature <temperature_default>` is used.

  :rotation:
    If the cell is filled with a universe, this element specifies the angles in
    degrees about the x, y, and z axes that the filled universe should be
    rotated. Should be given as three real numbers. For example, if you wanted
    to rotate the filled universe by 90 degrees about the z-axis, the cell
    element would look something like:

    .. code-block:: xml

        <cell fill="..." rotation="0 0 90" />

    The rotation applied is an intrinsic rotation whose Tait-Bryan angles are
    given as those specified about the x, y, and z axes respectively. That is to
    say, if the angles are :math:`(\phi, \theta, \psi)`, then the rotation
    matrix applied is :math:`R_z(\psi) R_y(\theta) R_x(\phi)` or

    .. math::

       \left [ \begin{array}{ccc} \cos\theta \cos\psi & -\cos\theta \sin\psi +
       \sin\phi \sin\theta \cos\psi & \sin\phi \sin\psi + \cos\phi \sin\theta
       \cos\psi \\ \cos\theta \sin\psi & \cos\phi \cos\psi + \sin\phi \sin\theta
       \sin\psi & -\sin\phi \cos\psi + \cos\phi \sin\theta \sin\psi \\
       -\sin\theta & \sin\phi \cos\theta & \cos\phi \cos\theta \end{array}
       \right ]

    *Default*: None

  :translation:
    If the cell is filled with a universe, this element specifies a vector that
    is used to translate (shift) the universe. Should be given as three real
    numbers.

    .. note:: Any translation operation is applied after a rotation, if also
              specified.

    *Default*: None


``<lattice>`` Element
---------------------

The ``<lattice>`` can be used to represent repeating structures (e.g. fuel pins
in an assembly) or other geometry which fits onto a rectilinear grid. Each cell
within the lattice is filled with a specified universe. A ``<lattice>`` accepts
the following attributes or sub-elements:

  :id:
    A unique integer that can be used to identify the lattice.

  :name:
    An optional string name to identify the lattice in summary output
    files. This string is limited to 52 characters for formatting purposes.

    *Default*: ""

  :dimension:
    Two or three integers representing the number of lattice cells in the x- and
    y- (and z-) directions, respectively.

    *Default*: None

  :lower_left:
    The coordinates of the lower-left corner of the lattice. If the lattice is
    two-dimensional, only the x- and y-coordinates are specified.

    *Default*: None

  :pitch:
    If the lattice is 3D, then three real numbers that express the distance
    between the centers of lattice cells in the x-, y-, and z- directions.  If
    the lattice is 2D, then omit the third value.

    *Default*: None

  :outer:
    The unique integer identifier of a universe that will be used to fill all
    space outside of the lattice.  The universe will be tiled repeatedly as if
    it were placed in a lattice of infinite size.  This element is optional.

    *Default*: An error will be raised if a particle leaves a lattice with no
    outer universe.

  :universes:
    A list of the universe numbers that fill each cell of the lattice.

    *Default*: None

Here is an example of a properly defined 2d rectangular lattice:

.. code-block:: xml

    <lattice id="10" dimension="3 3" outer="1">
        <lower_left> -1.5 -1.5 </lower_left>
        <pitch> 1.0 1.0 </pitch>
        <universes>
          2 2 2
          2 1 2
          2 2 2
        </universes>
    </lattice>

``<hex_lattice>`` Element
-------------------------

The ``<hex_lattice>`` can be used to represent repeating structures (e.g. fuel
pins in an assembly) or other geometry which naturally fits onto a hexagonal
grid or hexagonal prism grid. Each cell within the lattice is filled with a
specified universe. This lattice uses the "flat-topped hexagon" scheme where two
of the six edges are perpendicular to the y-axis.  A ``<hex_lattice>`` accepts
the following attributes or sub-elements:

  :id:
    A unique integer that can be used to identify the lattice.

  :name:
    An optional string name to identify the hex_lattice in summary output
    files. This string is limited to 52 characters for formatting purposes.

    *Default*: ""

  :n_rings:
    An integer representing the number of radial ring positions in the xy-plane.
    Note that this number includes the degenerate center ring which only has one
    element.

    *Default*: None

  :n_axial:
    An integer representing the number of positions along the z-axis.  This
    element is optional.

    *Default*: None

  :center:
    The coordinates of the center of the lattice. If the lattice does not have
    axial sections then only the x- and y-coordinates are specified.

    *Default*: None

  :pitch:
    If the lattice is 3D, then two real numbers that express the distance
    between the centers of lattice cells in the xy-plane and along the z-axis,
    respectively.  If the lattice is 2D, then omit the second value.

    *Default*: None

  :outer:
    The unique integer identifier of a universe that will be used to fill all
    space outside of the lattice.  The universe will be tiled repeatedly as if
    it were placed in a lattice of infinite size.  This element is optional.

    *Default*: An error will be raised if a particle leaves a lattice with no
    outer universe.

  :universes:
    A list of the universe numbers that fill each cell of the lattice.

    *Default*: None

Here is an example of a properly defined 2d hexagonal lattice:

.. code-block:: xml

    <hex_lattice id="10" n_rings="3" outer="1">
        <center> 0.0 0.0 </center>
        <pitch> 1.0 </pitch>
        <universes>
                  202
               202   202
            202   202   202
               202   202
            202   101   202
               202   202
            202   202   202
               202   202
                  202
        </universes>
    </hex_lattice>

.. _constructive solid geometry: http://en.wikipedia.org/wiki/Constructive_solid_geometry

.. _quadratic surfaces: http://en.wikipedia.org/wiki/Quadric

----------------------------------------
Materials Specification -- materials.xml
----------------------------------------

.. _cross_sections:

``<cross_sections>`` Element
----------------------------

The ``<cross_sections>`` element has no attributes and simply indicates the path
to an XML cross section listing file (usually named cross_sections.xml). If this
element is absent from the settings.xml file, the
:envvar:`OPENMC_CROSS_SECTIONS` environment variable will be used to find the
path to the XML cross section listing when in continuous-energy mode, and the
:envvar:`OPENMC_MG_CROSS_SECTIONS` environment variable will be used in
multi-group mode.

.. _multipole_library:

``<multipole_library>`` Element
-------------------------------

The ``<multipole_library>`` element indicates the directory containing a
windowed multipole library. If a windowed multipole library is available,
OpenMC can use it for on-the-fly Doppler-broadening of resolved resonance range
cross sections. If this element is absent from the settings.xml file, the
:envvar:`OPENMC_MULTIPOLE_LIBRARY` environment variable will be used.

  .. note:: The <temperature_multipole> element must also be set to "true" for
    windowed multipole functionality.

.. _material:

``<material>`` Element
----------------------

Each ``material`` element can have the following attributes or sub-elements:

  :id:
    A unique integer that can be used to identify the material.

  :name:
    An optional string name to identify the material in summary output
    files. This string is limited to 52 characters for formatting purposes.

    *Default*: ""

  :temperature:
    An element with no attributes which is used to set the default temperature
    of the material in Kelvin.

    *Default*: If a material default temperature is not given and a cell
    temperature is not specified, the :ref:`global default temperature
    <temperature_default>` is used.

  :density:
    An element with attributes/sub-elements called ``value`` and ``units``. The
    ``value`` attribute is the numeric value of the density while the ``units``
    can be "g/cm3", "kg/m3", "atom/b-cm", "atom/cm3", or "sum". The "sum" unit
    indicates that values appearing in ``ao`` or ``wo`` attributes for ``<nuclide>``
    and ``<element>`` sub-elements are to be interpreted as absolute nuclide/element
    densities in atom/b-cm or g/cm3, and the total density of the material is
    taken as the sum of all nuclides/elements. The "macro" unit is used with
    a ``macroscopic`` quantity to indicate that the density is already included
    in the library and thus not needed here.  However, if a value is provided
    for the ``value``, then this is treated as a number density multiplier on
    the macroscopic cross sections in the multi-group data.  This can be used,
    for example, when perturbing the density slightly.

    *Default*: None

    .. note:: A ``macroscopic`` quantity can not be used in conjunction with a
              ``nuclide``, ``element``, or ``sab`` quantity.

  :nuclide:
    An element with attributes/sub-elements called ``name``, and ``ao``
    or ``wo``. The ``name`` attribute is the name of the cross-section for a
    desired nuclide. Finally, the ``ao`` and ``wo`` attributes specify the atom or
    weight percent of that nuclide within the material, respectively. One
    example would be as follows:

    .. code-block:: xml

        <nuclide name="H1" ao="2.0" />
        <nuclide name="O16" ao="1.0" />

    .. note:: If one nuclide is specified in atom percent, all others must also
              be given in atom percent. The same applies for weight percentages.

    An optional attribute/sub-element for each nuclide is ``scattering``. This
    attribute may be set to "data" to use the scattering laws specified by the
    cross section library (default). Alternatively, when set to "iso-in-lab",
    the scattering laws are used to sample the outgoing energy but an
    isotropic-in-lab  distribution is used to sample the outgoing angle at each
    scattering interaction. The ``scattering`` attribute may be most useful
    when using OpenMC to compute multi-group cross-sections for deterministic
    transport codes and to quantify the effects of anisotropic scattering.

    *Default*: None

    .. note:: The ``scattering`` attribute/sub-element is not used in the
              multi-group :ref:`energy_mode`.

  :sab:
    Associates an S(a,b) table with the material. This element has one
    attribute/sub-element called ``name``. The ``name`` attribute
    is the name of the S(a,b) table that should be associated with the material.

    *Default*: None

    .. note:: This element is not used in the multi-group :ref:`energy_mode`.

  :macroscopic:
    The ``macroscopic`` element is similar to the ``nuclide`` element, but,
    recognizes that some multi-group libraries may be providing material
    specific macroscopic cross sections instead of always providing nuclide
    specific data like in the continuous-energy case.  To that end, the
    macroscopic element has one attribute/sub-element called ``name``.
    The ``name`` attribute is the name of the cross-section for a
    desired nuclide. One example would be as follows:

    .. code-block:: xml

        <macroscopic name="UO2" />

    .. note:: This element is only used in the multi-group :ref:`energy_mode`.

    *Default*: None

------------------------------------
Tallies Specification -- tallies.xml
------------------------------------

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

``<tally>`` Element
-------------------

The ``<tally>`` element accepts the following sub-elements:

  :name:
    An optional string name to identify the tally in summary output
    files. This string is limited to 52 characters for formatting purposes.

    *Default*: ""

  :filter:
    Specify a filter that restricts contributions to the tally to particles
    within certain regions of phase space. This element and its
    attributes/sub-elements are described below.

    .. note::
        You may specify zero, one, or multiple filters to apply to the tally. To
        specify multiple filters, you must use multiple ``<filter>`` elements.

    The ``filter`` element has the following attributes/sub-elements:

      :type:
        The type of the filter. Accepted options are "cell", "cellborn",
        "material", "universe", "energy", "energyout", "mesh", "distribcell",
        and "delayedgroup".

      :bins:
        For each filter type, the corresponding ``bins`` entry is given as
        follows:

        :cell:
          A list of cells in which the tally should be accumulated.

        :cellborn:
          This filter allows the tally to be scored to only when particles were
          originally born in a specified cell.

        :surface:
          A list of surfaces for which the tally should be accumulated.

        :material:
          A list of materials for which the tally should be accumulated.

        :universe:
          A list of universes for which the tally should be accumulated.

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
          The ``id`` of a structured mesh to be tallied over.

        :distribcell:
          The single cell which should be tallied uniquely for all instances.

          .. note::
              The distribcell filter will take a single cell ID and will tally
              each unique occurrence of that cell separately. This filter will
              not accept more than one cell ID. It is not recommended to combine
              this filter with a cell or mesh filter.

        :delayedgroup:
          A list of delayed neutron precursor groups for which the tally should
          be accumulated. For instance, to tally to all 6 delayed groups in the
          ENDF/B-VII.1 library the filter is specified as:

          .. code-block:: xml

              <filter type="delayedgroup" bins="1 2 3 4 5 6" />

          .. note:: This filter type is not used in the multi-group :ref:`energy_mode`.

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
        |                      |fission. This score type is not used in the        |
        |                      |multi-group :ref:`energy_mode`.                    |
        +----------------------+---------------------------------------------------+
        |prompt-nu-fission     |Total production of prompt neutrons due to         |
        |                      |fission. This score type is not used in the        |
        |                      |multi-group :ref:`energy_mode`.                    |
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
        |                      |This score type is not used in the                 |
        |                      |multi-group :ref:`energy_mode`.                    |
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
        |                      |This score type is not used in the                 |
        |                      |multi-group :ref:`energy_mode`.                    |
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

.. _usersguide_plotting:

--------------------------------------------
Geometry Plotting Specification -- plots.xml
--------------------------------------------

Basic plotting capabilities are available in OpenMC by creating a plots.xml
file and subsequently running with the command-line flag ``-plot``. The root
element of the plots.xml is simply ``<plots>`` and any number output plots can
be defined with ``<plot>`` sub-elements.  Two plot types are currently
implemented in openMC:

* ``slice``  2D pixel plot along one of the major axes. Produces a PPM image
  file.
* ``voxel``  3D voxel data dump. Produces a binary file containing voxel xyz
  position and cell or material id.


``<plot>`` Element
------------------

Each plot is specified by a combination of the following attributes or
sub-elements:

  :id:
    The unique ``id`` of the plot.

    *Default*: None - Required entry

  :filename:
    Filename for the output plot file.

    *Default*: "plot"

  :color:
    Keyword for plot coloring.  This can only be either ``cell`` or ``mat``,
    which colors regions by cells and materials, respectively. For voxel plots,
    this determines which id (cell or material) is associated with each
    position.

    *Default*: ``cell``

  :level:
    Universe depth to plot at (optional).  This parameter controls how many
    universe levels deep to pull cell and material ids from when setting plot
    colors.  If a given location does not have as many levels as specified,
    colors will be taken from the lowest level at that location. For example, if
    ``level`` is set to zero colors will be taken from top-level (universe zero)
    cells only.  However, if ``level`` is set to 1 colors will be taken from
    cells in universes that fill top-level fill-cells, and from top-level cells
    that contain materials.

    *Default*: Whatever the deepest universe is in the model

  :origin:
    Specifies the (x,y,z) coordinate of the center of the plot.  Should be three
    floats separated by spaces.

    *Default*: None - Required entry

  :width:
    Specifies the width of the plot along each of the basis directions.  Should
    be two or three floats separated by spaces for 2D plots and 3D plots,
    respectively.

    *Default*: None - Required entry

  :type:
    Keyword for type of plot to be produced. Currently only "slice" and "voxel"
    plots are implemented. The "slice" plot type creates 2D pixel maps saved in
    the PPM file format. PPM files can be displayed in most viewers (e.g. the
    default Gnome viewer, IrfanView, etc.).  The "voxel" plot type produces a
    binary datafile containing voxel grid positioning and the cell or material
    (specified by the ``color`` tag) at the center of each voxel. These
    datafiles can be processed into 3D SILO files using the
    ``openmc-voxel-to-silovtk`` utility provided with the OpenMC source, and
    subsequently viewed with a 3D viewer such as VISIT or Paraview. See the
    :ref:`io_voxel` for information about the datafile structure.

    .. note:: Since the PPM format is saved without any kind of compression,
              the resulting file sizes can be quite large.  Saving the image in
              the PNG format can often times reduce the file size by orders of
              magnitude without any loss of image quality. Likewise,
              high-resolution voxel files produced by OpenMC can be quite large,
              but the equivalent SILO files will be significantly smaller.

    *Default*: "slice"

``<plot>`` elements of ``type`` "slice" and "voxel" must contain the ``pixels``
attribute or sub-element:

  :pixels:
    Specifies the number of pixels or voxels to be used along each of the basis
    directions for "slice" and "voxel" plots, respectively. Should be two or
    three integers separated by spaces.

    .. warning:: The ``pixels`` input determines the output file size.  For the
                 PPM format, 10 million pixels will result in a file just under
                 30 MB in size. A 10 million voxel binary file will be around
                 40 MB.

    .. warning:: If the aspect ratio defined in ``pixels`` does not match the
                 aspect ratio defined in ``width`` the plot may appear stretched
                 or squeezed.

    .. warning:: Geometry features along a basis direction smaller than
                 ``width``/``pixels`` along that basis direction may not appear
                 in the plot.

    *Default*: None - Required entry for "slice" and "voxel" plots

``<plot>`` elements of ``type`` "slice" can also contain the following
attributes or sub-elements.  These are not used in "voxel" plots:

  :basis:
    Keyword specifying the plane of the plot for "slice" type plots.  Can be
    one of: "xy", "xz", "yz".

    *Default*: "xy"

  :background:
    Specifies the RGB color of the regions where no OpenMC cell can be found.
    Should be three integers separated by spaces.

    *Default*: 0 0 0 (black)

  :col_spec:
    Any number of this optional tag may be included in each ``<plot>`` element,
    which can override the default random colors for cells or materials. Each
    ``col_spec`` element must contain ``id`` and ``rgb`` sub-elements.

    :id:
      Specifies the cell or material unique id for the color specification.

    :rgb:
      Specifies the custom color for the cell or material. Should be 3 integers
      separated by spaces.

    As an example, if your plot is colored by material and you want material 23
    to be blue, the corresponding ``col_spec`` element would look like:

    .. code-block:: xml

        <col_spec id="23" rgb="0 0 255" />

    *Default*: None

  :mask:
    The special ``mask`` sub-element allows for the selective plotting of *only*
    user-specified cells or materials. Only one ``mask`` element is allowed per
    ``plot`` element, and it must contain as attributes or sub-elements a
    background masking color and a list of cells or materials to plot:

    :components:
      List of unique ``id`` numbers of the cells or materials to plot. Should be
      any number of integers separated by spaces.

    :background:
      Color to apply to all cells or materials not in the ``components`` list of
      cells or materials to plot. This overrides any ``col_spec`` color
      specifications.

    *Default*: None

  :meshlines:
    The ``meshlines`` sub-element allows for plotting the boundaries of a
    regular mesh on top of a plot. Only one ``meshlines`` element is allowed per
    ``plot`` element, and it must contain as attributes or sub-elements a mesh
    type and a linewidth.  Optionally, a color may be specified for the overlay:

    :meshtype:
      The type of the mesh to be plotted. Valid options are "tally", "entropy",
      "ufs", and "cmfd".  If plotting "tally" meshes, the id of the mesh to plot
      must be specified with the ``id`` sub-element.

    :id:
      A single integer id number for the mesh specified on ``tallies.xml`` that
      should be plotted. This element is only required for ``meshtype="tally"``.

    :linewidth:
      A single integer number of pixels of linewidth to specify for the mesh
      boundaries. Specifying this as 0 indicates that lines will be 1 pixel
      thick, specifying 1 indicates 3 pixels thick, specifying 2 indicates
      5 pixels thick, etc.

    :color:
      Specifies the custom color for the meshlines boundaries. Should be 3
      integers separated by whitespace.  This element is optional.

      *Default*: 0 0 0 (black)

    *Default*: None

.. _usersguide_cmfd:

------------------------------
CMFD Specification -- cmfd.xml
------------------------------

Coarse mesh finite difference acceleration method has been implemented in
OpenMC. Currently, it allows users to accelerate fission source convergence
during inactive neutron batches. To run CMFD, the ``<run_cmfd>`` element in
``settings.xml`` should be set to "true".

``<begin>`` Element
-------------------

The ``<begin>`` element controls what batch CMFD calculations should begin.

  *Default*: 1

``<dhat_reset>`` Element
------------------------

The ``<dhat_reset>`` element controls whether :math:`\widehat{D}` nonlinear
CMFD parameters should be reset to zero before solving CMFD eigenproblem.
It can be turned on with "true" and off with "false".

  *Default*: false

``<display>`` Element
---------------------

The ``<display>`` element sets one additional CMFD output column. Options are:

* "balance" - prints the RMS [%] of the resdiual from the neutron balance
  equation on CMFD tallies.
* "dominance" - prints the estimated dominance ratio from the CMFD iterations.
  **This will only work for power iteration eigensolver**.
* "entropy" - prints the *entropy* of the CMFD predicted fission source.
  **Can only be used if OpenMC entropy is active as well**.
* "source" - prints the RMS [%] between the OpenMC fission source and CMFD
  fission source.

  *Default*: balance

``<downscatter>`` Element
-------------------------

The ``<downscatter>`` element controls whether an effective downscatter cross
section should be used when using 2-group CMFD. It can be turned on with "true"
and off with "false".

  *Default*: false

``<feedback>`` Element
----------------------

The ``<feedback>`` element controls whether or not the CMFD diffusion result is
used to adjust the weight of fission source neutrons on the next OpenMC batch.
It can be turned on with "true" and off with "false".

  *Default*: false

``<gauss_seidel_tolerance>`` Element
------------------------------------

The ``<gauss_seidel_tolerance>`` element specifies two parameters. The first is
the absolute inner tolerance for Gauss-Seidel iterations when performing CMFD
and the second is the relative inner tolerance for Gauss-Seidel iterations
for CMFD calculations.

  *Default*: 1.e-10 1.e-5

``<ktol>`` Element
--------------------

The ``<ktol>`` element specifies the tolerance on the eigenvalue when performing
CMFD power iteration.

  *Default*: 1.e-8

``<mesh>`` Element
------------------

The CMFD mesh is a structured Cartesian mesh. This element has the following
attributes/sub-elements:

  :lower_left:
    The lower-left corner of the structured mesh. If only two coordinates are
    given, it is assumed that the mesh is an x-y mesh.

  :upper_right:
    The upper-right corner of the structrued mesh. If only two coordinates are
    given, it is assumed that the mesh is an x-y mesh.

  :dimension:
    The number of mesh cells in each direction.

  :width:
    The width of mesh cells in each direction.

  :energy:
    Energy bins [in eV], listed in ascending order (e.g. 0.0 0.625 20.0e6)
    for CMFD tallies and acceleration. If no energy bins are listed, OpenMC
    automatically assumes a one energy group calculation over the entire
    energy range.

  :albedo:
    Surface ratio of incoming to outgoing partial currents on global boundary
    conditions. They are listed in the following order: -x +x -y +y -z +z.

    *Default*: 1.0 1.0 1.0 1.0 1.0 1.0

  :map:
    An optional acceleration map can be specified to overlay on the coarse
    mesh spatial grid. If this option is used, a ``1`` is used for a
    non-accelerated region and a ``2`` is used for an accelerated region.
    For a simple 4x4 coarse mesh with a 2x2 fuel lattice surrounded by
    reflector, the map is:

      ``1 1 1 1``

      ``1 2 2 1``

      ``1 2 2 1``

      ``1 1 1 1``

    Therefore a 2x2 system of equations is solved rather than a 4x4. This
    is extremely important to use in reflectors as neutrons will not
    contribute to any tallies far away from fission source neutron regions.
    A ``2`` must be used to identify any fission source region.

    .. note:: Only two of the following three sub-elements are needed:
              ``lower_left``, ``upper_right`` and ``width``. Any combination
              of two of these will yield the third.

``<norm>`` Element
------------------

The ``<norm>`` element is used to normalize the CMFD fission source distribution
to a particular value. For example, if a fission source is calculated for a
17 x 17 lattice of pins, the fission source may be normalized to the number of
fission source regions, in this case 289. This is useful when visualizing this
distribution as the average peaking factor will be unity. This parameter will
not impact the calculation.

  *Default*: 1.0

``<power_monitor>`` Element
---------------------------

The ``<power_monitor>`` element is used to view the convergence of power
iteration. This option can be turned on with "true" and turned off with "false".

  *Default*: false

``<run_adjoint>`` Element
-------------------------

The ``<run_adjoint>`` element can be turned on with "true" to have an adjoint
calculation be performed on the last batch when CMFD is active.

  *Default*: false

``<shift>`` Element
--------------------

The ``<shift>`` element specifies an optional Wielandt shift parameter for
accelerating power iterations. It is by default very large so the impact of the
shift is effectively zero.

  *Default*: 1e6

``<spectral>`` Element
----------------------

The ``<spectral>`` element specifies an optional spectral radius that can be set to
accelerate the convergence of Gauss-Seidel iterations during CMFD power iteration
solve.

  *Default*: 0.0

``<stol>`` Element
------------------

The ``<stol>`` element specifies the tolerance on the fission source when performing
CMFD power iteration.

  *Default*: 1.e-8

``<tally_reset>`` Element
-------------------------

The ``<tally_reset>`` element contains a list of batch numbers in which CMFD tallies
should be reset.

  *Default*: None

``<write_matrices>`` Element
----------------------------

The ``<write_matrices>`` element is used to write the sparse matrices created
when solving CMFD equations. This option can be turned on with "true" and off
with "false".

  *Default*: false

------------------------------------
ERSN-OpenMC Graphical User Interface
------------------------------------

A third-party Java-based user-friendly graphical user interface for creating XML
input files called ERSN-OpenMC_ is developed and maintained by members of the
Radiation and Nuclear Systems Group at the Faculty of Sciences Tetouan, Morocco.
The GUI also allows one to automatically download prerequisites for installing and
running OpenMC.

.. _ERSN-OpenMC: https://github.com/EL-Bakkali-Jaafar/ERSN-OpenMC
