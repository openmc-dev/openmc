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

------------
Run Strategy
------------

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

.. _usersguide_batches:

Number of Batches
-----------------

In general, the stochastic uncertainty in your simulation results is directly
related to how many total active particles are simulated (the product of the
number of active batches, number of generations per batch, and number of
particles). At a minimum, you should use enough active batches so that the
central limit theorem is satisfied (about 30). Otherwise, reducing the overall
uncertainty in your simulation by a factor of 2 will require using 4 times as
many batches (since the standard deviation decreases as :math:`1/\sqrt{N}`).

Number of Inactive Batches
--------------------------

For :math:`k` eigenvalue simulations, the source distribution is not known a
priori. Thus, a "guess" of the source distribution is made and then iterated on,
with the source evolving closer to the true distribution at each iteration. Once
the source distribution has converged, it is then safe to start accumulating
tallies. Consequently, a preset number of inactive batches are run before the
active batches (where tallies are turned on) begin. The number of inactive
batches necessary to reach a converged source depends on the spatial extent of
the problem, its dominance ratio, what boundary conditions are used, and many
other factors. For small problems, using 50--100 inactive batches is likely
sufficient. For larger models, many hundreds of inactive batches may be
necessary. Users are recommended to use the :ref:`Shannon entropy
<usersguide_entropy>` diagnostic as a way of determining how many inactive
batches are necessary.

Specifying the initial source used for the very first batch is described in
:ref:`below <usersguide_source>`. Although the initial source is arbitrary in
the sense that any source will eventually converge to the correct distribution,
using a source guess that is closer to the actual converged source distribution
will translate into needing fewer inactive batches (and hence less simulation
time).

For fixed source simulations, the source distribution is known exactly, so no
inactive batches are needed. In this case the :attr:`Settings.inactive`
attribute can be omitted since it defaults to zero.

Number of Generations per Batch
-------------------------------

The standard deviation of tally results is calculated assuming that all
realizations (batches) are independent. However, in a :math:`k` eigenvalue
calculation, the source sites for each batch are produced from fissions in the
preceding batch, resulting in a correlation between successive batches. This
correlation can result in an underprediction of the variance. That is, the
variance reported is actually less than the true variance. To mitigate this
effect, OpenMC allows you to group together multiple fission generations into a
single batch for statistical purposes, rather than having each fission
generation be a separate batch, which is the default behavior.

Number of Particles per Generation
----------------------------------

There are several considerations for choosing the number of particles per
generation. As discussed in :ref:`usersguide_batches`, the total number of
active particles will determine the level of stochastic uncertainty in
simulation results, so using a higher number of particles will result in less
uncertainty. For parallel simulations that use OpenMP and/or MPI, the number of
particles per generation should be large enough to ensure good load balancing
between threads. For example, if you are running on a single processor with 32
cores, each core should have at least 100 particles or so (i.e., at least 3,200
particles per generation should be used). Using a larger number of particles per
generation can also help reduce the cost of synchronization and communication
between batches. For :math:`k` eigenvalue calculations, experts recommend_ at
least 10,000 particles per generation to avoid any bias in the estimate of
:math:`k` eigenvalue or tallies.

.. _recommend: https://permalink.lanl.gov/object/tr?what=info:lanl-repo/lareport/LA-UR-09-03136

.. _usersguide_source:

-----------------------------
External Source Distributions
-----------------------------

External source distributions can be specified through the
:attr:`Settings.source` attribute. If you have a single external source, you can
create an instance of any of the subclasses of :class:`openmc.SourceBase`
(:class:`openmc.IndependentSource`, :class:`openmc.FileSource`,
:class:`openmc.CompiledSource`) and use it to set the :attr:`Settings.source`
attribute. If you have multiple external sources with varying source strengths,
:attr:`Settings.source` should be set to a list of :class:`openmc.SourceBase`
objects.

The :class:`openmc.IndependentSource` class is the primary class for defining
source distributions and has four main attributes that one can set:
:attr:`IndependentSource.space`, which defines the spatial distribution,
:attr:`IndependentSource.angle`, which defines the angular distribution,
:attr:`IndependentSource.energy`, which defines the energy distribution, and
:attr:`IndependentSource.time`, which defines the time distribution.

The spatial distribution can be set equal to a sub-class of
:class:`openmc.stats.Spatial`; common choices are :class:`openmc.stats.Point` or
:class:`openmc.stats.Box`. To independently specify distributions in the
:math:`x`, :math:`y`, and :math:`z` coordinates, you can use
:class:`openmc.stats.CartesianIndependent`. To independently specify
distributions using spherical or cylindrical coordinates, you can use
:class:`openmc.stats.SphericalIndependent` or
:class:`openmc.stats.CylindricalIndependent`, respectively. Meshes can also be
used to represent spatial distributions with :class:`openmc.stats.MeshSpatial`
by specifying a mesh and source strengths for each mesh element. It is also
possible to define a "cloud" of source points, each with a different relative
probability, using :class:`openmc.stats.PointCloud`.

The angular distribution can be set equal to a sub-class of
:class:`openmc.stats.UnitSphere` such as :class:`openmc.stats.Isotropic`,
:class:`openmc.stats.Monodirectional`, or :class:`openmc.stats.PolarAzimuthal`.
By default, if no angular distribution is specified, an isotropic angular
distribution is used. As an example of a non-trivial angular distribution, the
following code would create a conical distribution with an aperture of 30
degrees pointed in the positive x direction::

  from math import pi, cos
  aperture = 30.0
  mu = openmc.stats.Uniform(cos(aperture/2), 1.0)
  phi = openmc.stats.Uniform(0.0, 2*pi)
  angle = openmc.stats.PolarAzimuthal(mu, phi, reference_uvw=(1., 0., 0.))

The energy distribution can be set equal to any univariate probability
distribution. This could be a probability mass function
(:class:`openmc.stats.Discrete`), a Watt fission spectrum
(:class:`openmc.stats.Watt`), or a tabular distribution
(:class:`openmc.stats.Tabular`). By default, if no energy distribution is
specified, a Watt fission spectrum with :math:`a` = 0.988 MeV and :math:`b` =
2.249 MeV :sup:`-1` is used.

The time distribution can be set equal to any univariate probability
distribution. This could be a probability mass function
(:class:`openmc.stats.Discrete`), a uniform distribution
(:class:`openmc.stats.Uniform`), or a tabular distribution
(:class:`openmc.stats.Tabular`). By default, if no time distribution is
specified, particles are started at :math:`t=0`.

As an example, to create an isotropic, 10 MeV monoenergetic source uniformly
distributed over a cube centered at the origin with an edge length of 10 cm, and
emitting a pulse of particles from 0 to 10 µs, one
would run::

  source = openmc.IndependentSource()
  source.space = openmc.stats.Box((-5, -5, -5), (5, 5, 5))
  source.angle = openmc.stats.Isotropic()
  source.energy = openmc.stats.Discrete([10.0e6], [1.0])
  source.time = openmc.stats.Uniform(0, 1e-6)
  settings.source = source

All subclasses of :class:`openmc.SourceBase` have a :attr:`SourceBase.strength`
attribute that indicates the relative strength of a source distribution if
multiple are used. For example, to create two sources, one that should be
sampled 70% of the time and another that should be sampled 30% of the time::

  src1 = openmc.IndependentSource()
  src1.strength = 0.7
  ...

  src2 = openmc.IndependentSource()
  src2.strength = 0.3
  ...

  settings.source = [src1, src2]

Finally, the :attr:`IndependentSource.particle` attribute can be used to
indicate the source should be composed of particles other than neutrons. For
example, the following would generate a photon source::

  source = openmc.IndependentSource()
  source.particle = 'photon'
  ...

  settings.source = source

For a full list of all classes related to statistical distributions, see
:ref:`pythonapi_stats`.

File-based Sources
------------------

OpenMC can use a pregenerated HDF5 source file through the
:class:`openmc.FileSource` class::

  settings.source = openmc.FileSource('source.h5')

Statepoint and source files are generated automatically when a simulation is run
and can be used as the starting source in a new simulation. Alternatively, a
source file can be manually generated with the :func:`openmc.write_source_file`
function. This is particularly useful for coupling OpenMC with another program
that generates a source to be used in OpenMC.

Surface Sources
+++++++++++++++

A source file based on particles that cross one or more surfaces can be
generated during a simulation using the :attr:`Settings.surf_source_write`
attribute::

  settings.surf_source_write = {
      'surfaces_ids': [1, 2, 3],
      'max_particles': 10000
  }

In this example, at most 10,000 source particles are stored when particles cross
surfaces with IDs of 1, 2, or 3. If no surface IDs are declared, particles
crossing any surface of the model will be banked::

  settings.surf_source_write = {'max_particles': 10000}

A cell ID can also be used to bank particles that are crossing any surface of
a cell that particles are either coming from or going to::

  settings.surf_source_write = {'cell': 1, 'max_particles': 10000}

In this example, particles that are crossing any surface that bounds cell 1 will
be banked excluding any surface that does not use a 'transmission' or 'vacuum'
boundary condition.

.. note:: Surfaces with boundary conditions that are not "transmission" or "vacuum"
          are not eligible to store any particles when using ``cell``, ``cellfrom``
          or ``cellto`` attributes. It is recommended to use surface IDs instead.

Surface IDs can be used in combination with a cell ID::

  settings.surf_source_write = {
      'cell': 1,
      'surfaces_ids': [1, 2, 3],
      'max_particles': 10000
  }

In that case, only particles that are crossing the declared surfaces coming from
cell 1 or going to cell 1 will be banked. To account specifically for particles
leaving or entering a given cell, ``cellfrom`` and ``cellto`` are also available
to respectively account for particles coming from a cell::

  settings.surf_source_write = {
      'cellfrom': 1,
      'max_particles': 10000
  }

or particles going to a cell::

  settings.surf_source_write = {
      'cellto': 1,
      'max_particles': 10000
  }

.. note:: The ``cell``, ``cellfrom`` and ``cellto`` attributes cannot be
          used simultaneously.

To generate more than one surface source files when the maximum number of stored
particles is reached, ``max_source_files`` is available. The surface source bank
will be cleared in simulation memory each time a surface source file is written.
As an example, to write a maximum of three surface source files:::

  settings.surf_source_write = {
      'surfaces_ids': [1, 2, 3],
      'max_particles': 10000,
      'max_source_files': 3
  }

.. _compiled_source:

Compiled Sources
----------------

It is often the case that one may wish to simulate a complex source distribution
that is not possible to represent with the classes described above. For these
situations, it is possible to define a complex source class containing an externally
defined source function that is loaded at runtime. A simple example source is shown
below.

.. code-block:: c++

  #include <memory> // for unique_ptr

  #include "openmc/random_lcg.h"
  #include "openmc/source.h"
  #include "openmc/particle.h"

  class CompiledSource : public openmc::Source
  {
    openmc::SourceSite sample(uint64_t* seed) const
    {
      openmc::SourceSite particle;
      // weight
      particle.particle = openmc::ParticleType::neutron;
      particle.wgt = 1.0;
      // position
      double angle = 2.0 * M_PI * openmc::prn(seed);
      double radius = 3.0;
      particle.r.x = radius * std::cos(angle);
      particle.r.y = radius * std::sin(angle);
      particle.r.z = 0.0;
      // angle
      particle.u = {1.0, 0.0, 0.0};
      particle.E = 14.08e6;
      particle.delayed_group = 0;
      return particle;
    }
  };

  extern "C" std::unique_ptr<CompiledSource> openmc_create_source(std::string parameters)
  {
    return std::make_unique<CompiledSource>();
  }

The above source creates monodirectional 14.08 MeV neutrons that are distributed
in a ring with a 3 cm radius. This routine is not particularly complex, but
should serve as an example upon which to build more complicated sources.

  .. note:: The source class must inherit from ``openmc::Source`` and
            implement a ``sample()`` function.

  .. note:: The ``openmc_create_source()`` function signature must be declared
            ``extern "C"``.

  .. note:: You should only use the ``openmc::prn()`` random number generator.

In order to build your external source, you will need to link it against the
OpenMC shared library. This can be done by writing a CMakeLists.txt file:

.. code-block:: cmake

   cmake_minimum_required(VERSION 3.3 FATAL_ERROR)
   project(openmc_sources CXX)
   add_library(source SHARED source_ring.cpp)
   find_package(OpenMC REQUIRED HINTS <path to openmc>)
   target_link_libraries(source OpenMC::libopenmc)

After running ``cmake`` and ``make``, you will have a libsource.so (or .dylib)
file in your build directory. You can then use this as an external source during
an OpenMC run by passing the path of the shared library to the
:class:`openmc.CompiledSource` class, which is then set as the
:attr:`Settings.source` attribute::

  settings.source = openmc.CompiledSource('libsource.so')

.. _parameterized_compiled_source:

Parameterized Compiled Sources
------------------------------

Some compiled sources may have values (parameters) that can be changed between
runs. This is supported by using the ``openmc_create_source()`` function to pass
parameters to the source class when it is created:

.. code-block:: c++

  #include <memory> // for unique_ptr

  #include "openmc/source.h"
  #include "openmc/particle.h"

  class CompiledSource : public openmc::Source {
  public:
    CompiledSource(double energy) : energy_{energy} { }

    // Samples from an instance of this class.
    openmc::SourceSite sample(uint64_t* seed) const
    {
      openmc::SourceSite particle;
      // weight
      particle.particle = openmc::ParticleType::neutron;
      particle.wgt = 1.0;
      // position
      particle.r.x = 0.0;
      particle.r.y = 0.0;
      particle.r.z = 0.0;
      // angle
      particle.u = {1.0, 0.0, 0.0};
      particle.E = this->energy_;
      particle.delayed_group = 0;

      return particle;
    }

  private:
    double energy_;
  };

  extern "C" std::unique_ptr<CompiledSource> openmc_create_source(std::string parameter) {
    double energy = std::stod(parameter);
    return std::make_unique<CompiledSource>(energy);
  }

When creating an instance of the :class:`openmc.CompiledSource` class, you will
need to pass both the path of the shared library as well as the parameters as a
string, which gets passed down to the ``openmc_create_source()`` function::

  settings.source = openmc.CompiledSource('libsource.so', '3.5e6')

.. _usersguide_source_constraints:

Source Constraints
------------------

All source classes in OpenMC have the ability to apply a set of "constraints"
that limit which sampled source sites are actually used for transport. The most
common use case is to sample source sites over some simple spatial distribution
(e.g., uniform over a box) and then only accept those that appear in a given
cell or material. This can be done with a domain constraint, which can be
specified as follows::

  source_cell = openmc.Cell(...)
  ...

  spatial_dist = openmc.stats.Box((-10., -10., -10.), (10., 10., 10.))
  source = openmc.IndependentSource(
      space=spatial_dist,
      constraints={'domains': [source_cell]}
  )

For k-eigenvalue problems, a convenient constraint is available that limits
source sites to those sampled in a fissionable material::

  source = openmc.IndependentSource(
      space=spatial_dist, constraints={'fissionable': True}
  )

Constraints can also be placed on a range of energies or times::

  # Only use source sites between 500 keV and 1 MeV and with times under 1 sec
  source = openmc.FileSource(
      'source.h5',
      constraints={'energy_bounds': [500.0e3, 1.0e6], 'time_bounds': [0.0, 1.0]}
  )

Normally, when a source site is rejected, a new one will be resampled until one
is found that meets the constraints. However, the rejection strategy can be
changed so that a rejected site will just not be simulated by specifying::

  source = openmc.IndependentSource(
      space=spatial_dist,
      constraints={'domains': [cell], 'rejection_strategy': 'kill'}
  )

In this case, the actual number of particles simulated may be less than what you
specified in :attr:`Settings.particles`.

.. _usersguide_entropy:

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

Particle Track Files
--------------------

OpenMC can generate a particle track file that contains track information
(position, direction, energy, time, weight, cell ID, and material ID) for each
state along a particle's history. There are two ways to indicate which particles
and/or how many particles should have their tracks written. First, you can
identify specific source particles by their batch, generation, and particle ID
numbers::

  settings.tracks = [
    (1, 1, 50),
    (2, 1, 30),
    (5, 1, 75)
  ]

In this example, track information would be written for the 50th particle in the
1st generation of batch 1, the 30th particle in the first generation of batch 2,
and the 75th particle in the 1st generation of batch 5. Unless you are using
more than one generation per batch (see :ref:`usersguide_particles`), the
generation number should be 1. Alternatively, you can run OpenMC in a mode where
track information is written for *all* particles, up to a user-specified limit::

  openmc.run(tracks=True)

In this case, you can control the maximum number of source particles for which
tracks will be written as follows::

  settings.max_tracks = 1000

Particle track information is written to the ``tracks.h5`` file, which can be
analyzed using the :class:`~openmc.Tracks` class::

  >>> tracks = openmc.Tracks('tracks.h5')
  >>> tracks
  [<Track (1, 1, 50): 151 particles>,
   <Track (2, 1, 30): 191 particles>,
   <Track (5, 1, 75): 81 particles>]

Each :class:`~openmc.Track` object stores a list of track information for every
primary/secondary particle. In the above example, the first source particle
produced 150 secondary particles for a total of 151 particles. Information for
each primary/secondary particle can be accessed using the
:attr:`~openmc.Track.particle_tracks` attribute::

  >>> first_track = tracks[0]
  >>> first_track.particle_tracks
  [<ParticleTrack: neutron, 120 states>,
   <ParticleTrack: photon, 6 states>,
   <ParticleTrack: electron, 2 states>,
   <ParticleTrack: electron, 2 states>,
   <ParticleTrack: electron, 2 states>,
   ...
   <ParticleTrack: electron, 2 states>,
   <ParticleTrack: electron, 2 states>]
  >>> photon = first_track.particle_tracks[1]

The :class:`~openmc.ParticleTrack` class is a named tuple indicating the
particle type and then a NumPy array of the "states". The states array is a
compound type with a field for each physical quantity (position, direction,
energy, time, weight, cell ID, and material ID). For example, to get the
position for the above particle track::

  >>> photon.states['r']
  array([(-11.92987939, -12.28467295, 0.67837495),
         (-11.95213726, -12.2682    , 0.68783964),
         (-12.2682    , -12.03428339, 0.82223855),
         (-12.5913778 , -11.79510096, 0.95966298),
         (-12.6622572 , -11.74264344, 0.98980293),
         (-12.6907775 , -11.7215357 , 1.00193058)],
        dtype=[('x', '<f8'), ('y', '<f8'), ('z', '<f8')])

The full list of fields is as follows:

  :r: Position (each direction in [cm])
  :u: Direction
  :E: Energy in [eV]
  :time: Time in [s]
  :wgt: Weight
  :cell_id: Cell ID
  :cell_instance: Cell instance
  :material_id: Material ID

Both the :class:`~openmc.Tracks` and :class:`~openmc.Track` classes have a
``filter`` method that allows you to get a subset of tracks that meet a given
criteria. For example, to get all tracks that involved a photon::

  >>> tracks.filter(particle='photon')
  [<Track (1, 1, 50): 151 particles>,
   <Track (2, 1, 30): 191 particles>,
   <Track (5, 1, 75): 81 particles>]

The :meth:`openmc.Tracks.filter` method returns a new :class:`~openmc.Tracks`
instance, whereas the :meth:`openmc.Track.filter` method returns a new
:class:`~openmc.Track` instance.

.. note:: If you are using an MPI-enabled install of OpenMC and run a simulation
          with more than one process, a separate track file will be written for
          each MPI process with the filename ``tracks_p#.h5`` where # is the
          rank of the corresponding process. Multiple track files can be
          combined with the :ref:`scripts_track_combine` script:

          .. code-block:: sh

            openmc-track-combine tracks_p*.h5 --out tracks.h5

-----------------------
Restarting a Simulation
-----------------------

OpenMC can be run in a mode where it reads in a statepoint file and continues a
simulation from the ending point of the statepoint file. A restart simulation
can be performed by passing the path to the statepoint file to the OpenMC
executable:

.. code-block:: sh

    openmc -r statepoint.100.h5

From the Python API, the `restart_file` argument provides the same behavior:

.. code-block:: python

    openmc.run(restart_file='statepoint.100.h5')

or if using the :class:`~openmc.Model` class:

.. code-block:: python

    model.run(restart_file='statepoint.100.h5')

The restart simulation will execute until the number of batches specified in the
:class:`~openmc.Settings` object on a model (or in the :ref:`settings XML file
<io_settings>`) is satisfied. Note that if the number of batches in the
statepoint file is the same as that specified in the settings object (i.e., if
the inputs were not modified before the restart run), no particles will be
transported and OpenMC will exit immediately.

.. note:: A statepoint file must match the input model to be successfully used in a restart simulation.
