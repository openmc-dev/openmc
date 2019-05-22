import numpy as np

from .angle_energy import AngleEnergy
from .correlated import CorrelatedAngleEnergy


class CoherentElasticAE(AngleEnergy):
    r"""Differential cross section for coherent elastic scattering

    The differential cross section for coherent elastic scattering from a
    powdered crystalline material may be represented as:

    .. math::
        \frac{d^2\sigma}{dE'd\Omega} (E\rightarrow E',\mu,T) = \frac{1}{E} \sum
        \limits_{i=1}^{E_i < E} s_i(T) \delta(\mu - \mu_i) \delta (E - E')
        /(2\pi)

    where :math:`E_i` are the energies of the Bragg edges in [eV], :math:`s_i(T)`
    is the structure factor in [eV-b] at the moderator temperature :math:`T`
    in [K], and :math:`\mu_i = 1 - 2E_i/E`.

    Parameters
    ----------
    coherent_xs : openmc.data.CoherentElastic
        Coherent elastic scattering cross section

    Attributes
    ----------
    coherent_xs : openmc.data.CoherentElastic
        Coherent elastic scattering cross section

    """
    def __init__(self, coherent_xs):
        self.coherent_xs = coherent_xs

    def to_hdf5(self, group):
        """Write coherent elastic distribution to an HDF5 group

        Parameters
        ----------
        group : h5py.Group
            HDF5 group to write to

        """
        group.attrs['type'] = np.string_('coherent_elastic')
        group['coherent_xs'] = group.parent['xs']


class IncoherentElasticAE(AngleEnergy):
    r"""Differential cross section for incoherent elastic scattering

    The differential cross section for incoherent elastic scattering may be
    represented as:

    .. math::
        \frac{d^2\sigma}{dE'd\Omega} (E\rightarrow E',\mu,T) = \frac{\sigma_b}
        {4\pi} e^{-2EW'(T)(1-\mu)} \delta(E - E')

    where :math:`\sigma_b` is the characteristic cross section in [b] and
    :math:`W'(T)` is the Debye-Waller integral divided by the atomic mass in
    [eV\ :math:`^{-1}`].

    Parameters
    ----------
    debye_waller : float
        Debye-Waller integral in [eV\ :math:`^{-1}`]

    Attributes
    ----------
    debye_waller : float
        Debye-Waller integral in [eV\ :math:`^{-1}`]

    """
    def __init__(self, debye_waller):
        self.debye_waller = debye_waller

    def to_hdf5(self, group):
        """Write incoherent elastic distribution to an HDF5 group

        Parameters
        ----------
        group : h5py.Group
            HDF5 group to write to

        """
        group.attrs['type'] = np.string_('incoherent_elastic')
        group.create_dataset('debye_waller', data=self.debye_waller)

    @classmethod
    def from_hdf5(cls, group):
        """Generate incoherent elastic distribution from HDF5 data

        Parameters
        ----------
        group : h5py.Group
            HDF5 group to read from

        Returns
        -------
        openmc.data.IncoherentElasticAE
            Incoherent elastic distribution

        """
        return cls(group['debye_waller'])


class IncoherentElasticAEDiscrete(AngleEnergy):
    """Discrete angle representation of incoherent elastic scattering

    Parameters
    ----------
    mu_out : numpy.ndarray
        Equi-probable discrete angles at each incoming energy

    """
    def __init__(self, mu_out):
        self.mu_out = mu_out

    def to_hdf5(self, group):
        """Write discrete incoherent elastic distribution to an HDF5 group

        Parameters
        ----------
        group : h5py.Group
            HDF5 group to write to

        """
        group.attrs['type'] = np.string_('incoherent_elastic_discrete')
        group.create_dataset('mu_out', data=self.mu_out)

    @classmethod
    def from_hdf5(cls, group):
        """Generate discrete incoherent elastic distribution from HDF5 data

        Parameters
        ----------
        group : h5py.Group
            HDF5 group to read from

        Returns
        -------
        openmc.data.IncoherentElasticAEDiscrete
            Discrete incoherent elastic distribution

        """
        return cls(group['mu_out'][()])


class IncoherentInelasticAEDiscrete(AngleEnergy):
    """Discrete angle representation of incoherent inelastic scattering

    Parameters
    ----------
    energy_out : numpy.ndarray
        Outgoing energies for each incoming energy
    mu_out : numpy.ndarray
        Discrete angles for each incoming/outgoing energy
    skewed : bool
        Whether discrete angles are equi-probable or have a skewed distribution

    Attributes
    ----------
    energy_out : numpy.ndarray
        Outgoing energies for each incoming energy
    mu_out : numpy.ndarray
        Discrete angles for each incoming/outgoing energy
    skewed : bool
        Whether discrete angles are equi-probable or have a skewed distribution

    """
    def __init__(self, energy_out, mu_out, skewed=False):
        self.energy_out = energy_out
        self.mu_out = mu_out
        self.skewed = skewed

    def to_hdf5(self, group):
        """Write discrete incoherent inelastic distribution to an HDF5 group

        Parameters
        ----------
        group : h5py.Group
            HDF5 group to write to

        """
        group.attrs['type'] = np.string_('incoherent_inelastic_discrete')
        group.create_dataset('energy_out', data=self.energy_out)
        group.create_dataset('mu_out', data=self.mu_out)
        group.create_dataset('skewed', data=self.skewed)

    @classmethod
    def from_hdf5(cls, group):
        """Generate discrete incoherent inelastic distribution from HDF5 data

        Parameters
        ----------
        group : h5py.Group
            HDF5 group to read from

        Returns
        -------
        openmc.data.IncoherentInelasticAEDiscrete
            Discrete incoherent inelastic distribution

        """
        energy_out = group['energy_out'][()]
        mu_out = group['mu_out'][()]
        skewed = bool(group['skewed'])
        return cls(energy_out, mu_out, skewed)


class IncoherentInelasticAE(CorrelatedAngleEnergy):
    _name = 'incoherent_inelastic'
