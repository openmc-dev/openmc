.. _variance_reduction:

==================
Variance Reduction
==================

Global variance reduction in OpenMC is accomplished by weight windowing
or source biasing techniques, the latter of which additionally provides a 
local variance reduction capability. OpenMC is capable of generating weight 
windows using either the MAGIC or FW-CADIS methods. Both techniques will 
produce a ``weight_windows.h5`` file that can be loaded and used later on. In 
this section, we first break down the steps required to generate and apply 
weight windows, then describe how source biasing may be applied.

.. _ww_generator:

------------------------------------
Generating Weight Windows with MAGIC
------------------------------------

As discussed in the :ref:`methods section <methods_variance_reduction>`, MAGIC
is an iterative method that uses flux tally information from a Monte Carlo
simulation to produce weight windows for a user-defined mesh. While generating
the weight windows, OpenMC is capable of applying the weight windows generated
from a previous batch while processing the next batch, allowing for progressive
improvement in the weight window quality across iterations.

The typical way of generating weight windows is to define a mesh and then add an
:class:`openmc.WeightWindowGenerator` object to an :attr:`openmc.Settings`
instance, as follows::

    # Define weight window spatial mesh
    ww_mesh = openmc.RegularMesh()
    ww_mesh.dimension = (10, 10, 10)
    ww_mesh.lower_left = (0.0, 0.0, 0.0)
    ww_mesh.upper_right = (100.0, 100.0, 100.0)

    # Create weight window object and adjust parameters
    wwg = openmc.WeightWindowGenerator(
        method='magic',
        mesh=ww_mesh,
        max_realizations=settings.batches
    )

    # Add generator to Settings instance
    settings.weight_window_generators = wwg

Notably, the :attr:`max_realizations` attribute is adjusted to the number of
batches, such that all iterations are used to refine the weight window
parameters.

With the :class:`~openmc.WeightWindowGenerator` instance added to the
:attr:`~openmc.Settings`, the rest of the problem can be defined as normal. When
running, note that the second iteration and beyond may be several orders of
magnitude slower than the first. As the weight windows are applied in each
iteration, particles may be agressively split, resulting in a large number of
secondary (split) particles being generated per initial source particle. This is
not necessarily a bad thing, as the split particles are much more efficient at
exploring low flux regions of phase space as compared to initial particles.
Thus, even though the reported "particles/second" metric of OpenMC may be much
lower when generating (or just applying) weight windows as compared to analog
MC, it typically leads to an overall improvement in the figure of merit
accounting for the reduction in the variance.

.. warning::
    The number of particles per batch may need to be adjusted downward
    significantly to result in reasonable runtimes when weight windows are being
    generated or used.

At the end of the simulation, a ``weight_windows.h5`` file will be saved to disk
for later use. Loading it in another subsequent simulation will be discussed in
the "Using Weight Windows" section below.

------------------------------------------------------
Generating Weight Windows with FW-CADIS and Random Ray
------------------------------------------------------

Weight window generation with FW-CADIS and random ray in OpenMC uses the same
exact strategy as with MAGIC. An :class:`openmc.WeightWindowGenerator` object is
added to the :attr:`openmc.Settings` object, and a ``weight_windows.h5`` will be
generated at the end of the simulation. The only difference is that the code
must be run in random ray mode. A full description of how to enable and setup
random ray mode can be found in the :ref:`Random Ray User Guide <random_ray>`.

.. note::
    It is a long term goal for OpenMC to be able to generate FW-CADIS weight
    windows with only a few tweaks to an existing continuous energy Monte Carlo
    input deck. However, at the present time, the workflow requires several
    steps to generate multigroup cross section data and to configure the random
    ray solver. A high level overview of the current workflow for generation of
    weight windows with FW-CADIS using random ray is given below.

1. Begin by making a deepy copy of your continuous energy Python model and then
   convert the copy to be multigroup and use the random ray transport solver.
   The conversion process can largely be automated as described in more detail
   in the :ref:`random ray quick start guide <quick_start>`, summarized below::

    # Define continuous energy model
    ce_model = openmc.pwr_pin_cell() # example, replace with your model

    # Make a copy to convert to multigroup and random ray
    model = copy.deepcopy(ce_model)

    # Convert model to multigroup (will auto-generate MGXS library if needed)
    model.convert_to_multigroup()

    # Convert model to random ray and initialize random ray parameters
    # to reasonable defaults based on the specifics of the geometry
    model.convert_to_random_ray()

    # (Optional) Overlay source region decomposition mesh to improve fidelity of the
    # random ray solver. Adjust 'n' for fidelity vs runtime.
    n = 10
    mesh = openmc.RegularMesh()
    mesh.dimension = (n, n, n)
    mesh.lower_left = model.geometry.bounding_box.lower_left
    mesh.upper_right = model.geometry.bounding_box.upper_right
    model.settings.random_ray['source_region_meshes'] = [(mesh, [model.geometry.root_universe])]

    # (Optional) Improve fidelity of the random ray solver by enabling linear sources
    model.settings.random_ray['source_shape'] = 'linear'

    # (Optional) Increase the number of rays/batch, to reduce uncertainty
    model.settings.particles = 500

   If you need to improve the fidelity of the MGXS library, there is more
   information on generating multigroup cross sections via OpenMC in the
   :ref:`random ray MGXS guide <mgxs_gen>`.

2. Add in a :class:`~openmc.WeightWindowGenerator` in a similar manner as for
   MAGIC generation with Monte Carlo and set the :attr:`method` attribute set to
   ``"fw_cadis"``::

       # Create weight window object and adjust parameters, using the same mesh
       # we used for source region decomposition
       wwg = openmc.WeightWindowGenerator(
           method='fw_cadis',
           mesh=mesh,
           max_realizations=settings.batches
       )

       # Add generator to openmc.settings object
       settings.weight_window_generators = wwg

.. warning::
    If using FW-CADIS weight window generation, ensure that the selected weight
    window mesh does not subdivide any source regions in the problem. This can
    be ensured by using the same mesh for both source region subdivision (i.e.,
    assigning to ``model.settings.random_ray['source_region_meshes']``) and for
    weight window generation.

3. When running your multigroup random ray input deck, OpenMC will automatically
   run a forward solve followed by an adjoint solve, with a
   ``weight_windows.h5`` file generated at the end. The ``weight_windows.h5``
   file will contain FW-CADIS generated weight windows. This file can be used in
   identical manner as one generated with MAGIC, as described below.

--------------------
Using Weight Windows
--------------------

To use a ``weight_windows.h5`` weight window file with OpenMC's Monte Carlo
solver, the Python input just needs to load the h5 file::

    settings.weight_window_checkpoints = {'collision': True, 'surface': True}
    settings.survival_biasing = False
    settings.weight_windows = openmc.WeightWindowsList.from_hdf5('weight_windows.h5')
    settings.weight_windows_on = True

The :class:`~openmc.WeightWindowGenerator` instance is not needed to load an
existing ``weight_windows.h5`` file. Inclusion of a
:class:`~openmc.WeightWindowGenerator` instance will cause OpenMC to generate
*new* weight windows and thus overwrite the existing ``weight_windows.h5`` file.
Weight window mesh information is embedded into the weight window file, so the
mesh does not need to be redefined. Monte Carlo solves that load a weight window
file as above will utilize weight windows to reduce the variance of the
simulation.

.. _source_biasing:

--------------
Source Biasing
--------------

In fixed source problems in OpenMC, source biasing provides a means to reduce 
the variance on global or localized responses, depending on the biasing scheme. 
In either case, the premise of the method is to sample particle source sites 
from a biasing distribution that directs a larger fraction of the simulated 
histories towards phase space regions of interest than would be found there 
under analog sampling. In order to preserve an unbiased estimate of the tally 
mean, the weight of these particles is adjusted according to the probability of 
sampling the source site with analog sampling, divided by the probability 
assigned by the bias distribution. While the assignment of statistical weights 
is outlined in the :ref:`methods section <methods_source_biasing>`, this 
section demonstrates the implementation of source biasing to problems in 
OpenMC.

Source biasing in OpenMC is accomplished by applying a distribution to the 
:attr:`bias` attribute of one or more of the univariate or independent 
multivariate distributions which make up an :class:`~openmc.IndependentSource` 
instance as follows::

    # First instantiate the bias distribution:
    bias_dist = openmc.stats.PowerLaw(a=0,b=3,n=3)

    # Construct a new distribution with the bias applied:
    biased_dist = openmc.stats.PowerLaw(a=0,b=3,n=2,bias=bias_dist)

    # The bias attribute can also be set on an existing "analog" distribution:
    sphere_dist = openmc.stats.spherical_uniform(r_outer=3) 
    sphere_dist.r.bias = bias_dist

Univariate distributions (:class:`openmc.stats.Discrete`, 
:class:`openmc.stats.Uniform`, :class:`openmc.stats.PowerLaw`, 
:class:`openmc.stats.Maxwell`, :class:`openmc.stats.Watt`, 
:class:`openmc.stats.Normal`, :class:`openmc.stats.Tabular`, and 
:class:`openmc.stats.Mixture` instances combining these) may be sampled via the 
Python API, returning the sample(s) along with the associated weight(s)::

    sample_vec, wgt_vec = biased_dist.sample(n_samples=100,seed=1)

Here, if the distribution is unbiased, the weight of each sample will be unity. 
Finally, :class:`~openmc.IndependentSource` instances can be constructed with 
biased distributions and added to an :class:`openmc.Settings` instance as 
usual::

    source = openmc.IndependentSource()
    source.space = sphere_dist # Biased sampling will occur for spatial coordinate
    source.angle = openmc.stats.Isotropic()
    source.energy = openmc.stats.Discrete([10.0e6], [1.0])
    source.time = openmc.stats.Uniform(0, 1e-6)
    settings.source = source

During the simulation, source sites are then sampled using the bias 
distributions where available, and given starting statistical weights 
corresponding to the cumulative product of the weights assigned by each 
distribution in the source object. Hence multiple source variables (e.g. 
direction and energy) may be biased and the resulting source sites will have 
their weights adjusted accordingly.

.. note::
    Combining source biasing with weight windows can be a powerful variance 
    reduction technique if each is constructed appropriately for the response
    of interest. For example, if a source biasing scheme is devised for 
    variance reduction of a specific localized response, the user may be able 
    to specify their own weight window structure that results in more efficient  
    transport than if weight windows were generated by either of OpenMC's 
    automatic weight window generators, which are intended for global variance 
    reduction. Consistency must also be enforced between the source particle 
    weights and the window into which they are born, so that unnecessary 
    splitting and roulette do not occur at the start of each history.

Biasing distributions which could result in degenerate weight mappings is not 
recommended; this is most commonly seen when biasing the :math:`\phi`-coordinate 
of spherical or cylindrical independent multivariate distributions. In such 
cases degenerate behavior will be observed at the pole about which :math:`\phi` 
is measured, with all values of :math:`\phi` (hence many possible statistical 
weights) mapping to the same point for :math:`r=0`` or :math:`\mu=0`, and 
large weight gradients in the vicinity. In most cases requiring a spherical 
independent source, it would be preferable to reorient the reference vector of 
the distribution such that biasing could be applied to the 
:math:`\mu`-coordinate instead.

When biasing a distribution, care should also be taken to ensure that both the 
unbiased and the biasing distribution share a common support--that is, every 
region of phase space mapped to a nonzero probability density by the unbiased 
distribution should likewise map to nonzero probability under the biased 
distribution, and vice-versa. In OpenMC, this places restrictions on the set of 
compatible distributions which may be used to bias sampling of each 
distribution type. The following table summarizes the method for each 
distribution in OpenMC which permits biasing.

.. table:: **Biasable probability density functions (PDFs) in OpenMC**

    +----------------------------------------------+---------------------------+
    |Discrete Univariate PDFs                      | Biasing Method            |
    +==============================================+===========================+
    | :class:`openmc.stats.Discrete`               |Apply a second, unbiased   |
    |                                         | :class:`openmc.stats.Discrete` |
    |                                              |sharing the same :attr:`x` |
    |                                              |vector to the :attr:`bias` |
    |                                              |attribute                  |
    +----------------------------------------------+---------------------------+

    +----------------------------------------------+---------------------------+
    |Continuous Univariate PDFs                    | Biasing Method            |
    +==============================================+===========================+
    | :class:`openmc.stats.Uniform`,               |Apply a second, unbiased   |
    | :class:`openmc.stats.PowerLaw`,              |continuous univariate PDF  |
    | :class:`openmc.stats.Maxwell`,               |to the :attr:`bias`        |
    | :class:`openmc.stats.Watt`,                  |attribute                  |
    | :class:`openmc.stats.Normal`,                |                           |
    | :class:`openmc.stats.Tabular`                |                           |
    +----------------------------------------------+---------------------------+

    +----------------------------------------------+---------------------------+
    |Discrete Multivariate PDFs                    | Biasing Method            |
    +==============================================+===========================+
    | :class:`openmc.stats.PointCloud`,            |Apply a vector of the new  |
    | :class:`openmc.stats.MeshSpatial`            |relative probabilities of  |
    |                                              |each point or mesh element |
    |                                              |under biased sampling to   |
    |                                              |the :attr:`bias` attribute |
    +----------------------------------------------+---------------------------+

    +----------------------------------------------+---------------------------+
    |Continuous Multivariate PDFs                  | Biasing Method            |
    +==============================================+===========================+
    | :class:`openmc.stats.CartesianIndependent`,  |Construct from biased      |
    | :class:`openmc.stats.CylindricalIndependent`,|univariate distributions   |
    | :class:`openmc.stats.SphericalIndependent`,  |for :attr:`x`, :attr:`y`,  | 
    | :class:`openmc.stats.PolarAzimuthal`         | :attr:`z`, etc.           |
    +----------------------------------------------+---------------------------+
    | :class:`openmc.stats.Isotropic`              |Apply an unbiased          |
    |                                   | :class:`openmc.stats.PolarAzimuthal` |
    |                                              | to the :attr:`bias`       |
    |                                              | attribute                 |
    +----------------------------------------------+---------------------------+

.. note::
    The :class:`openmc.stats.Mixture` class may be constructed from multiple 
    biased univariate distributions, but biasing the probabilities used to 
    select between these distributions is not currently supported.
