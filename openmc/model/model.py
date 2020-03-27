from collections.abc import Iterable
from pathlib import Path
import time

import openmc
from openmc.checkvalue import check_type, check_value


class Model:
    """Model container.

    This class can be used to store instances of :class:`openmc.Geometry`,
    :class:`openmc.Materials`, :class:`openmc.Settings`,
    :class:`openmc.Tallies`, :class:`openmc.Plots`, and :class:`openmc.CMFD`,
    thus making a complete model. The :meth:`Model.export_to_xml` method will
    export XML files for all attributes that have been set. If the
    :meth:`Model.materials` attribute is not set, it will attempt to create a
    ``materials.xml`` file based on all materials appearing in the geometry.

    Parameters
    ----------
    geometry : openmc.Geometry, optional
        Geometry information
    materials : openmc.Materials, optional
        Materials information
    settings : openmc.Settings, optional
        Settings information
    tallies : openmc.Tallies, optional
        Tallies information
    plots : openmc.Plots, optional
        Plot information

    Attributes
    ----------
    geometry : openmc.Geometry
        Geometry information
    materials : openmc.Materials
        Materials information
    settings : openmc.Settings
        Settings information
    tallies : openmc.Tallies
        Tallies information
    plots : openmc.Plots
        Plot information

    """

    def __init__(self, geometry=None, materials=None, settings=None,
                 tallies=None, plots=None):
        self.geometry = openmc.Geometry()
        self.materials = openmc.Materials()
        self.settings = openmc.Settings()
        self.tallies = openmc.Tallies()
        self.plots = openmc.Plots()

        if geometry is not None:
            self.geometry = geometry
        if materials is not None:
            self.materials = materials
        if settings is not None:
            self.settings = settings
        if tallies is not None:
            self.tallies = tallies
        if plots is not None:
            self.plots = plots

    @property
    def geometry(self):
        return self._geometry

    @property
    def materials(self):
        return self._materials

    @property
    def settings(self):
        return self._settings

    @property
    def tallies(self):
        return self._tallies

    @property
    def plots(self):
        return self._plots

    @geometry.setter
    def geometry(self, geometry):
        check_type('geometry', geometry, openmc.Geometry)
        self._geometry = geometry

    @materials.setter
    def materials(self, materials):
        check_type('materials', materials, Iterable, openmc.Material)
        if isinstance(materials, openmc.Materials):
            self._materials = materials
        else:
            del self._materials[:]
            for mat in materials:
                self._materials.append(mat)

    @settings.setter
    def settings(self, settings):
        check_type('settings', settings, openmc.Settings)
        self._settings = settings

    @tallies.setter
    def tallies(self, tallies):
        check_type('tallies', tallies, Iterable, openmc.Tally)
        if isinstance(tallies, openmc.Tallies):
            self._tallies = tallies
        else:
            del self._tallies[:]
            for tally in tallies:
                self._tallies.append(tally)

    @plots.setter
    def plots(self, plots):
        check_type('plots', plots, Iterable, openmc.Plot)
        if isinstance(plots, openmc.Plots):
            self._plots = plots
        else:
            del self._plots[:]
            for plot in plots:
                self._plots.append(plot)

    def deplete(self, timesteps, chain_file=None, method='cecm',
                fission_q=None, **kwargs):
        """Deplete model using specified timesteps/power

        Parameters
        ----------
        timesteps : iterable of float
            Array of timesteps in units of [s]. Note that values are not
            cumulative.
        chain_file : str, optional
            Path to the depletion chain XML file.  Defaults to the chain
            found under the ``depletion_chain`` in the
            :envvar:`OPENMC_CROSS_SECTIONS` environment variable if it exists.
        method : str
             Integration method used for depletion (e.g., 'cecm', 'predictor')
        fission_q : dict, optional
            Dictionary of nuclides and their fission Q values [eV].
            If not given, values will be pulled from the ``chain_file``.
        **kwargs
            Keyword arguments passed to integration function (e.g.,
            :func:`openmc.deplete.integrator.cecm`)

        """
        # Import the depletion module.  This is done here rather than the module
        # header to delay importing openmc.lib (through openmc.deplete) which
        # can be tough to install properly.
        import openmc.deplete as dep

        # Create OpenMC transport operator
        op = dep.Operator(
            self.geometry, self.settings, chain_file,
            fission_q=fission_q,
        )

        # Perform depletion
        check_value('method', method, ('cecm', 'predictor', 'cf4', 'epc_rk4',
                                       'si_celi', 'si_leqi', 'celi', 'leqi'))
        getattr(dep.integrator, method)(op, timesteps, **kwargs)

    def export_to_xml(self, directory='.'):
        """Export model to XML files.

        Parameters
        ----------
        directory : str
            Directory to write XML files to. If it doesn't exist already, it
            will be created.

        """
        # Create directory if required
        d = Path(directory)
        if not d.is_dir():
            d.mkdir(parents=True)

        self.settings.export_to_xml(d)
        if not self.settings.dagmc:
            self.geometry.export_to_xml(d)

        # If a materials collection was specified, export it. Otherwise, look
        # for all materials in the geometry and use that to automatically build
        # a collection.
        if self.materials:
            self.materials.export_to_xml(d)
        else:
            materials = openmc.Materials(self.geometry.get_all_materials()
                                         .values())
            materials.export_to_xml(d)

        if self.tallies:
            self.tallies.export_to_xml(d)
        if self.plots:
            self.plots.export_to_xml(d)

    def run(self, **kwargs):
        """Creates the XML files, runs OpenMC, and returns the path to the last
        statepoint file generated.

        Parameters
        ----------
        **kwargs
            Keyword arguments passed to :func:`openmc.run`

        Returns
        -------
        Path
            Path to the last statepoint written by this run
            (None if no statepoint was written)

        """

        self.export_to_xml()

        # Setting tstart here ensures we don't pick up any pre-existing statepoint
        # files in the output directory
        tstart = time.time()
        last_statepoint = None

        openmc.run(**kwargs)

        # Get output directory and return the last statepoint written by this run
        if self.settings.output and 'path' in self.settings.output:
            output_dir = Path(self.settings.output['path'])
        else:
            output_dir = Path.cwd()
        for sp in output_dir.glob('statepoint.*.h5'):
            mtime = sp.stat().st_mtime
            if mtime >= tstart:  # >= allows for poor clock resolution
                tstart = mtime
                last_statepoint = sp
        return last_statepoint
