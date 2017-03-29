import openmc
from openmc.checkvalue import check_type


class Model(object):
    """OpenMC model container for the openmc.Geometry, openmc.Materials,
    openmc.Settings, openmc.Tallies, openmc.CMFD objects, and openmc.Plot
    objects

    Parameters
    ----------
    geometry : openmc.Geometry
        Geometry information
    materials : openmc.Materials
        Materials information
    settings : openmc.Settings
        Settings information
    tallies : openmc.Tallies
        Tallies information, optional
    cmfd : openmc.CMFD
        CMFD information, optional
    plots : openmc.Plots
        Plot information, optional

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
    cmfd : openmc.CMFD
        CMFD information
    plots : openmc.Plots
        Plot information

    """

    def __init__(self, geometry, materials, settings, tallies=None, cmfd=None,
                 plots=None):
        self.geometry = geometry
        self.materials = materials
        self.settings = settings
        if tallies:
            self.tallies = tallies
        else:
            self._tallies = openmc.Tallies()
        if cmfd:
            self.cmfd = cmfd
        else:
            self._cmfd = None
        if plots:
            self.plots = plots
        else:
            self.plots = openmc.Plots()

        self.sp = None

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
    def cmfd(self):
        return self._cmfd

    @property
    def plots(self):
        return self._plots

    @geometry.setter
    def geometry(self, geometry):
        check_type('geometry', geometry, openmc.Geometry)
        self._geometry = geometry

    @materials.setter
    def materials(self, materials):
        check_type('materials', materials, openmc.Materials)
        self._materials = materials

    @settings.setter
    def settings(self, settings):
        check_type('settings', settings, openmc.Settings)
        self._settings = settings

    @tallies.setter
    def tallies(self, tallies):
        check_type('tallies', tallies, openmc.Tallies)
        self._tallies = tallies

    @cmfd.setter
    def cmfd(self, cmfd):
        check_type('cmfd', cmfd, openmc.CMFD)
        self._cmfd = cmfd

    @plots.setter
    def plots(self, plots):
        check_type('plots', plots, openmc.Plots)
        self._plots = plots

    def export_to_xml(self):
        """Export model settings to XML files.
        """

        self.geometry.export_to_xml()
        self.materials.export_to_xml()
        self.settings.export_to_xml()
        self.tallies.export_to_xml()
        if self.cmfd is not None:
            self.cmfd.export_to_xml()
        self.plots.export_to_xml()

    def run(self, **kwargs):
        """Creates the XML files, runs OpenMC, and loads the statepoint.

        Parameters
        ----------
        **kwargs
            All keyword arguments are passed to openmc.run

        Returns
        -------
        2-tuple of float
            k_combined from the statepoint

        """

        self.export_to_xml()

        return_code = openmc.run(**kwargs)

        assert (return_code == 0), "OpenMC did not execute successfully"

        statepoint_batches = self.settings.batches
        if self.settings.statepoint is not None:
            if 'batches' in self.settings.statepoint:
                statepoint_batches = self.settings.statepoint['batches'][-1]
        self.sp = \
            openmc.StatePoint('statepoint.{}.h5'.format(statepoint_batches))

        return self.sp.k_combined

    def close(self):
        """Close the statepoint and summary files
        """

        if self.sp is not None:
            self.sp._f.close()
            self.sp.summary._f.close()
