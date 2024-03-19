from pathlib import Path

import openmc


def test_from_endf_w774():
    """W74 HDF5 data."""
    filename = Path('endf-b-vii.1') / 'protons' / 'p-080_Hg_200.endf'
    openmc.data.IncidentChargedParticle.from_endf(filename)