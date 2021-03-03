.. _usersguide_volume:

==============================
Stochastic Volume Calculations
==============================

.. currentmodule:: openmc

OpenMC has a capability to stochastically determine volumes of cells, materials,
and universes. The method works by overlaying a bounding box, sampling points
from within the box, and seeing what fraction of points were found in a desired
domain. The benefit of doing this stochastically (as opposed to equally-spaced
points), is that it is possible to give reliable error estimates on each
stochastic quantity.

To specify that a volume calculation be run, you first need to create an
instance of :class:`openmc.VolumeCalculation`. The constructor takes a list of
cells, materials, or universes; the number of samples to be used; and the
lower-left and upper-right Cartesian coordinates of a bounding box that encloses
the specified domains::

  lower_left = (-0.62, -0.62, -50.)
  upper_right = (0.62, 0.62, 50.)
  vol_calc = openmc.VolumeCalculation([fuel, clad, moderator], 1000000,
                                      lower_left, upper_right)

For domains contained within regions that have simple definitions, OpenMC can
sometimes automatically determine a bounding box. In this case, the last two
arguments are not necessary. For example,

::

   sphere = openmc.Sphere(r=10.0)
   cell = openm.Cell(region=-sphere)
   vol_calc = openmc.VolumeCalculation([cell], 1000000)

Of course, the volumes that you *need* this capability for are often the ones
with complex definitions.

A threshold can be applied for the calculation's variance, standard deviation,
or relative error of volume estimates using :meth:`openmc.VolumeCalculation.set_trigger`::

    vol_calc.set_trigger(1e-05, 'std_dev')

If a threshold is provided, calculations will be performed iteratively using the
number of samples specified on the calculation until all volume uncertainties are below
the threshold value. If no threshold is provided, the calculation will run the number of
samples specified once and return the result.

Once you have one or more :class:`openmc.VolumeCalculation` objects created, you
can then assign then to :attr:`Settings.volume_calculations`::

  settings = openmc.Settings()
  settings.volume_calculations = [cell_vol_calc, mat_vol_calc]

To execute the volume calculations, one can either set :attr:`Settings.run_mode`
to 'volume' and run :func:`openmc.run`, or alternatively run
:func:`openmc.calculate_volumes` which doesn't require that
:attr:`Settings.run_mode` be set.

When your volume calculations have finished, you can load the results using the
:meth:`VolumeCalculation.load_results` method on an existing object. If you
don't have an existing :class:`VolumeCalculation` object, you can create one and
load results simultaneously using the :meth:`VolumeCalculation.from_hdf5` class
method::

  vol_calc = openmc.VolumeCalculation(...)
  ...
  openmc.calculate_volumes()
  vol_calc.load_results('volume_1.h5')

  # ..or we can create a new object
  vol_calc = openmc.VolumeCalculation.from_hdf5('volume_1.h5')

After the results are loaded, volume estimates will be stored in
:attr:`VolumeCalculation.volumes`. There is also a
:attr:`VolumeCalculation.atoms_dataframe` attribute that shows stochastic
estimates of the number of atoms of each type of nuclide within the specified
domains along with their uncertainties.
