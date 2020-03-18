.. _usersguide_settings:

==================
Execution Settings
==================

.. currentmodule:: openmc

Once you have created the materials and geometry for your simulation, the last
step to have a complete model is to specify execution settings through the
:class:`openmc.Settings` class. At a minimum, you need to specify a :ref:`source
distribution <usersguide_source>` and :ref:`how many particles to run
<usersguide_particles>`. Many other execution settings can be set using the
:class:`openmc.Settings` object, but they are generally optional.

.. _usersguide_run_modes:

---------
Run Modes
---------

The :attr:`Settings.run_mode` attribute controls what run mode is used when
:ref:`scripts_openmc` is executed. There are five different run modes that can
be specified:

'eigenvalue'
  Runs a :math:`k` eigenvalue simulation. See :ref:`methods_eigenvalue` for a
  full description of eigenvalue calculations. In this mode, the
  :attr:`Settings.source` specifies a starting source that is only used for the
  first fission generation.

'fixed source'
  Runs a fixed-source calculation with a specified external source, specified in
  the :attr:`Settings.source` attribute.

'volume'
  Runs a stochastic volume calculation.

'plot'
  Generates slice or voxel plots (see :ref:`usersguide_plots`).

'particle restart'
  Simulate a single source particle using a particle restart file.


So, for example, to specify that OpenMC should be run in fixed source mode, you
would need to instantiate a :class:`openmc.Settings` object and assign the
:attr:`Settings.run_mode` attribute::

   settings = openmc.Settings()
   settings.run_mode = 'fixed source'

If you don't specify a run mode, the default run mode is 'eigenvalue'.

.. _usersguide_particles:

-------------------
Number of Particles
-------------------

For a fixed source simulation, the total number of source particle histories
simulated is broken up into a number of *batches*, each corresponding to a
:ref:`realization <methods_tallies>` of the tally random variables. Thus, you
need to specify both the number of batches (:attr:`Settings.batches`) as well as
the number of particles per batch (:attr:`Settings.particles`).

For a :math:`k` eigenvalue simulation, particles are grouped into *fission
generations*, as described in :ref:`methods_eigenvalue`. Successive fission
generations can be combined into a batch for statistical purposes. By default, a
batch will consist of only a single fission generation, but this can be changed
with the :attr:`Settings.generations_per_batch` attribute. For problems with a
high dominance ratio, using multiple generations per batch can help reduce
underprediction of variance, thereby leading to more accurate confidence
intervals. Tallies should not be scored to until the source distribution
converges, as described in :ref:`method-successive-generations`, which may take
many generations. To specify the number of batches that should be discarded
before tallies begin to accumulate, use the :attr:`Settings.inactive` attribute.

The following example shows how one would simulate 10000 particles per
generation, using 10 generations per batch, 150 total batches, and discarding 5
batches. Thus, a total of 145 active batches (or 1450 generations) will be used
for accumulating tallies.

::

   settings.particles = 10000
   settings.generations_per_batch = 10
   settings.batches = 150
   settings.inactive = 5

.. _usersguide_source:

-----------------------------
External Source Distributions
-----------------------------

External source distributions can be specified through the
:attr:`Settings.source` attribute. If you have a single external source, you can
create an instance of :class:`openmc.Source` and use it to set the
:attr:`Settings.source` attribute. If you have multiple external sources with
varying source strengths, :attr:`Settings.source` should be set to a list of
:class:`openmc.Source` objects.

The :class:`openmc.Source` class has three main attributes that one can set:
:attr:`Source.space`, which defines the spatial distribution,
:attr:`Source.angle`, which defines the angular distribution, and
:attr:`Source.energy`, which defines the energy distribution.

The spatial distribution can be set equal to a sub-class of
:class:`openmc.stats.Spatial`; common choices are :class:`openmc.stats.Point` or
:class:`openmc.stats.Box`. To independently specify distributions in the
:math:`x`, :math:`y`, and :math:`z` coordinates, you can use
:class:`openmc.stats.CartesianIndependent`.

The angular distribution can be set equal to a sub-class of
:class:`openmc.stats.UnitSphere` such as :class:`openmc.stats.Isotropic`,
:class:`openmc.stats.Monodirectional`, or
:class:`openmc.stats.PolarAzimuthal`. By default, if no angular distribution is
specified, an isotropic angular distribution is used.

The energy distribution can be set equal to any univariate probability
distribution. This could be a probability mass function
(:class:`openmc.stats.Discrete`), a Watt fission spectrum
(:class:`openmc.stats.Watt`), or a tabular distribution
(:class:`openmc.stats.Tabular`). By default, if no energy distribution is
specified, a Watt fission spectrum with :math:`a` = 0.988 MeV and :math:`b` =
2.249 MeV :sup:`-1` is used.

As an example, to create an isotropic, 10 MeV monoenergetic source uniformly
distributed over a cube centered at the origin with an edge length of 10 cm, one
would run::

  source = openmc.Source()
  source.space = openmc.stats.Box((-5, -5, -5), (5, 5, 5))
  source.angle = openmc.stats.Isotropic()
  source.energy = openmc.stats.Discrete([10.0e6], [1.0])
  settings.source = source

The :class:`openmc.Source` class also has a :attr:`Source.strength` attribute
that indicates the relative strength of a source distribution if multiple are
used. For example, to create two sources, one that should be sampled 70% of the
time and another that should be sampled 30% of the time::

  src1 = openmc.Source()
  src1.strength = 0.7
  ...

  src2 = openmc.Source()
  src2.strength = 0.3
  ...

  settings.source = [src1, src2]

Finally, the :attr:`Source.particle` attribute can be used to indicate the
source should be composed of particles other than neutrons. For example, the
following would generate a photon source::

  source = openmc.Source()
  source.particle = 'photon'
  ...

  settings.source = source

For a full list of all classes related to statistical distributions, see
:ref:`pythonapi_stats`.

---------------
Shannon Entropy
---------------

To assess convergence of the source distribution, the scalar Shannon entropy
metric is often used in Monte Carlo codes. OpenMC also allows you to calculate
Shannon entropy at each generation over a specified mesh, created using the
:class:`openmc.RegularMesh` class. After instantiating a :class:`RegularMesh`,
you need to specify the lower-left coordinates of the mesh
(:attr:`RegularMesh.lower_left`), the number of mesh cells in each direction
(:attr:`RegularMesh.dimension`) and either the upper-right coordinates of the
mesh (:attr:`RegularMesh.upper_right`) or the width of each mesh cell
(:attr:`RegularMesh.width`). Once you have a mesh, simply assign it to the
:attr:`Settings.entropy_mesh` attribute.

::

   entropy_mesh = openmc.RegularMesh()
   entropy_mesh.lower_left = (-50, -50, -25)
   entropy_mesh.upper_right = (50, 50, 25)
   entropy_mesh.dimension = (8, 8, 8)

   settings.entropy_mesh = entropy_mesh

If you're unsure of what bounds to use for the entropy mesh, you can try getting
a bounding box for the entire geometry using the :attr:`Geometry.bounding_box`
property::

  geom = openmc.Geometry()
  ...
  m = openmc.RegularMesh()
  m.lower_left, m.upper_right = geom.bounding_box
  m.dimension = (8, 8, 8)

  settings.entropy_mesh = m

----------------
Photon Transport
----------------

In addition to neutrons, OpenMC is also capable of simulating the passage of
photons through matter. This allows the modeling of photon production from
neutrons as well as pure photon calculations. The
:attr:`Settings.photon_transport` attribute can be used to enable photon
transport::

  settings.photon_transport = True

The way in which OpenMC handles secondary charged particles can be specified
with the :attr:`Settings.electron_treatment` attribute. By default, the
:ref:`thick-target bremsstrahlung <ttb>` (TTB) approximation is used to generate
bremsstrahlung radiation emitted by electrons and positrons created in photon
interactions. To neglect secondary bremsstrahlung photons and instead deposit
all energy from electrons locally, the local energy deposition option can be
selected::

  settings.electron_treatment = 'led'

.. note::
   Some features related to photon transport are not currently implemented,
   including:

     * Tallying photon energy deposition.
     * Generating a photon source from a neutron calculation that can be used
       for a later fixed source photon calculation.
     * Photoneutron reactions.

--------------------------
Generation of Output Files
--------------------------

A number of attributes of the :class:`openmc.Settings` class can be used to
control what files are output and how often. First, there is the
:attr:`Settings.output` attribute which takes a dictionary having keys
'summary', 'tallies', and 'path'. The first two keys controls whether a
``summary.h5`` and ``tallies.out`` file are written, respectively (see
:ref:`result_files` for a description of those files). By default, output files
are written to the current working directory; this can be changed by setting the
'path' key. For example, if you want to disable the ``tallies.out`` file and
write the ``summary.h5`` to a directory called 'results', you'd specify the
:attr:`Settings.output` dictionary as::

   settings.output = {
       'tallies': False,
       'path': 'results'
   }

Generation of statepoint and source files is handled separately through the
:attr:`Settings.statepoint` and :attr:`Settings.sourcepoint` attributes. Both of
those attributes expect dictionaries and have a 'batches' key which indicates at
which batches statepoints and source files should be written. Note that by
default, the source is written as part of the statepoint file; this behavior can
be changed by the 'separate' and 'write' keys of the
:attr:`Settings.sourcepoint` dictionary, the first of which indicates whether
the source should be written to a separate file and the second of which
indicates whether the source should be written at all.

As an example, to write a statepoint file every five batches::

  settings.batches = n
  settings.statepoint = {'batches': range(5, n + 5, 5)}

.. _NIST ESTAR database: https://physics.nist.gov/PhysRefData/Star/Text/ESTAR.html
