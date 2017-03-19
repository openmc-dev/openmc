import openmc
from openmc.checkvalue import check_type


class Model(object):
    """OpenMC model container for the openmc.Geometry, openmc.Materials,
    openmc.Settings, openmc.Tallies, and openmc.CMFD objects

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

    """

    def __init__(self, geometry, materials, settings, tallies=None, cmfd=None):
        self.geometry = geometry
        self.materials = materials
        self.settings = settings
        if tallies:
            self.tallies = tallies
        else:
            self._tallies = None
        if cmfd:
            self.cmfd = cmfd
        else:
            self._cmfd = None

        self.sp = None

    @property
    def geometry(self):
        return self._geometry

    @geometry.setter
    def geometry(self, geometry):
        check_type('geometry', geometry, openmc.Geometry)
        self._geometry = geometry

    @property
    def materials(self):
        return self._materials

    @materials.setter
    def materials(self, materials):
        check_type('materials', materials, openmc.Materials)
        self._materials = materials

    @property
    def settings(self):
        return self._settings

    @settings.setter
    def settings(self, settings):
        check_type('settings', settings, openmc.Settings)
        self._settings = settings

    @property
    def tallies(self):
        return self._tallies

    @tallies.setter
    def tallies(self, tallies):
        check_type('tallies', tallies, openmc.Tallies)
        self._tallies = tallies

    @property
    def cmfd(self):
        return self._cmfd

    @cmfd.setter
    def cmfd(self, cmfd):
        check_type('cmfd', cmfd, openmc.CMFD)
        self._cmfd = cmfd

    def export_to_xml(self):
        """Export model settings to XML files.
        """

        self.geometry.export_to_xml()
        self.materials.export_to_xml()
        self.settings.export_to_xml()
        if self.tallies:
            self.tallies.export_to_xml()
        if self.cmfd:
            self.cmfd.export_to_xml()

    def execute(self, output=True):
        """Creates the XML files, runs OpenMC, and loads the statepoint.

        Parameters
        ----------
        output : bool
            Capture OpenMC output from standard out, defaults to True

        Returns
        -------
        Iterable of float
            k_combined from the statepoint

        """

        self.export_to_xml()

        openmc.run(output=output)

        statepoint_batches = self.settings.batches
        if self.settings.statepoint is not None:
            if 'batches' in self.settings.statepoint:
                statepoint_batches = self.settings.statepoint['batches'][-1]
        self.sp = openmc.StatePoint('statepoint.' + str(statepoint_batches) +
                                    '.h5')

        return self.sp.k_combined

    def close(self):
        """Close the statepoint and summary files
        """

        if self.sp is not None:
            self.sp._f.close()
            self.sp.summary._f.close()
