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
create an instance of :class:`openmc.Source` and use it to set the
:attr:`Settings.source` attribute. If you have multiple external sources with
varying source strengths, :attr:`Settings.source` should be set to a list of
:class:`openmc.Source` objects.

The :class:`openmc.Source` class has four main attributes that one can set:
:attr:`Source.space`, which defines the spatial distribution,
:attr:`Source.angle`, which defines the angular distribution,
:attr:`Source.energy`, which defines the energy distribution, and
:attr:`Source.time`, which defines the time distribution.

The spatial distribution can be set equal to a sub-class of
:class:`openmc.stats.Spatial`; common choices are :class:`openmc.stats.Point` or
:class:`openmc.stats.Box`. To independently specify distributions in the
:math:`x`, :math:`y`, and :math:`z` coordinates, you can use
:class:`openmc.stats.CartesianIndependent`. To independently specify
distributions using spherical or cylindrical coordinates, you can use
:class:`openmc.stats.SphericalIndependent` or
:class:`openmc.stats.CylindricalIndependent`, respectively.

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

  source = openmc.Source()
  source.space = openmc.stats.Box((-5, -5, -5), (5, 5, 5))
  source.angle = openmc.stats.Isotropic()
  source.energy = openmc.stats.Discrete([10.0e6], [1.0])
  source.time = openmc.stats.Uniform(0, 1e-6)
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

File-based Sources
------------------

OpenMC can use a pregenerated HDF5 source file by specifying the ``filename``
argument to :class:`openmc.Source`::

  settings.source = openmc.Source(filename='source.h5')

Statepoint and source files are generated automatically when a simulation is run
and can be used as the starting source in a new simulation. Alternatively, a
source file can be manually generated with the :func:`openmc.write_source_file`
function. This is particularly useful for coupling OpenMC with another program
that generates a source to be used in OpenMC.

A source file based on particles that cross one or more surfaces can be
generated during a simulation using the :attr:`Settings.surf_source_write`
attribute::

  settings.surf_source_write = {
      'surfaces_ids': [1, 2, 3],
      'max_particles': 10000
  }

In this example, at most 10,000 source particles are stored when particles cross
surfaces with IDs of 1, 2, or 3.

.. _custom_source:

Custom Sources
--------------

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

  class CustomSource : public openmc::Source
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

  extern "C" std::unique_ptr<CustomSource> openmc_create_source(std::string parameters)
  {
    return std::make_unique<CustomSource>();
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
file in your build directory. Setting the :attr:`openmc.Source.library`
attribute to the path of this shared library will indicate that it should be
used for sampling source particles at runtime.

.. _parameterized_custom_source:

Custom Parameterized Sources
----------------------------

Some custom sources may have values (parameters) that can be changed between
runs. This is supported by using the ``openmc_create_source()`` function to
pass parameters defined in the :attr:`openmc.Source.parameters` attribute to
the source class when it is created:

.. code-block:: c++

  #include <memory> // for unique_ptr

  #include "openmc/source.h"
  #include "openmc/particle.h"

  class CustomSource : public openmc::Source {
  public:
    CustomSource(double energy) : energy_{energy} { }

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

  extern "C" std::unique_ptr<CustomSource> openmc_create_source(std::string parameter) {
    double energy = std::stod(parameter);
    return std::make_unique<CustomSource>(energy);
  }

As with the basic custom source functionality, the custom source library
location must be provided in the :attr:`openmc.Source.library` attribute.

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
