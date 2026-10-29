"""Test that the extra neutrons of (n,2n) reactions keep the lifetime clock
that the IFP generation time is scored with."""

import math

import openmc

from tests.regression_tests import config

# Distance flown through void by the source neutrons [cm]
FLIGHT_PATH = 1.0e4
# Energy of the source neutrons [eV]
SOURCE_ENERGY = 1.4e7
# Speed of light [cm/s] and neutron rest mass energy [eV]
SPEED_OF_LIGHT = 2.99792458e10
NEUTRON_MASS = 9.3956542052e8


def n2n_clock_model():
    """Eigenvalue model in which every neutron of the first generation
    descends from a 14 MeV neutron that crosses 100 m of void before it
    reaches a beryllium and uranium target."""
    beryllium = openmc.Material()
    beryllium.add_nuclide('Be9', 1.0)
    beryllium.set_density('g/cm3', 1.85)
    uranium = openmc.Material()
    uranium.add_nuclide('U235', 1.0)
    uranium.set_density('g/cm3', 16.0)

    side = openmc.model.RectangularPrism(30.0, 30.0, boundary_type='vacuum')
    planes = [
        openmc.ZPlane(-1.0, boundary_type='vacuum'),
        openmc.ZPlane(FLIGHT_PATH),
        openmc.ZPlane(FLIGHT_PATH + 5.0),
        openmc.ZPlane(FLIGHT_PATH + 15.0, boundary_type='vacuum'),
    ]
    model = openmc.Model()
    model.geometry = openmc.Geometry([
        openmc.Cell(region=-side & +planes[0] & -planes[1]),
        openmc.Cell(region=-side & +planes[1] & -planes[2], fill=beryllium),
        openmc.Cell(region=-side & +planes[2] & -planes[3], fill=uranium),
    ])

    model.settings.particles = 5000
    model.settings.batches = 2
    model.settings.inactive = 1
    model.settings.ifp_n_generation = 1
    model.settings.source = openmc.IndependentSource(
        space=openmc.stats.Point(),
        angle=openmc.stats.Monodirectional((0.0, 0.0, 1.0)),
        energy=openmc.stats.Discrete([SOURCE_ENERGY], [1.0]),
        time=openmc.stats.Discrete([1.0], [1.0]),
    )

    tally = openmc.Tally()
    tally.scores = ['ifp-time-numerator', 'ifp-denominator']
    model.tallies = [tally]
    return model


def test_n2n_neutrons_keep_lifetime_clock(run_in_tmpdir):
    """Check the IFP mean lifetime against the flight time of the source
    neutrons.

    Fission neutrons start the next generation, and nothing is split, so a
    neutron of the first generation is either a source neutron or an extra
    neutron of an (n,2n) or other (n,xn) reaction with an integral yield,
    which continues the history of the neutron that produced it. Every
    source neutron flies 100 m through void before its first collision, and
    a neutron that leaves the target for the void cannot come back, since
    the void is convex. If the extra neutrons keep the lifetime clock of the
    neutron that produced them, every fission site of the first generation
    is made by a neutron whose clock is at least the flight time. The IFP
    tallies of the second generation score the lifetimes of those sites, so
    the ratio of the time numerator to the denominator, a mean over
    positive weights, is at least the flight time too, whatever the random
    numbers. The bound uses the non-relativistic speed, which exceeds the
    relativistic speed of the transport, so it is below the actual flight
    time. If the extra neutrons restarted the clock at the reaction, in the
    target, their fission sites would record some nanoseconds and the ratio
    would fall below the bound. The source neutrons start at t = 1 s, so a
    clock that copied the time instead of the lifetime would exceed the
    upper bound.
    """
    model = n2n_clock_model()
    kwargs = {'openmc_exec': config['exe'], 'event_based': config['event']}
    if config['mpi']:
        kwargs['mpi_args'] = [config['mpiexec'], '-n', config['mpi_np']]
    sp_path = model.run(**kwargs)

    tally = model.tallies[0]
    with openmc.StatePoint(sp_path) as sp:
        result = sp.tallies[tally.id]
        numerator = result.get_values(scores=['ifp-time-numerator']).item()
        denominator = result.get_values(scores=['ifp-denominator']).item()

    speed = SPEED_OF_LIGHT * math.sqrt(2.0 * SOURCE_ENERGY / NEUTRON_MASS)
    flight_time = FLIGHT_PATH / speed
    assert denominator > 0.0
    mean_lifetime = numerator / denominator
    assert mean_lifetime >= flight_time, (
        f'IFP mean lifetime ({mean_lifetime:.4e} s) is below the flight '
        f'time of the source neutrons ({flight_time:.4e} s)'
    )
    assert mean_lifetime <= 1.0e-3, (
        f'IFP mean lifetime ({mean_lifetime:.4e} s) exceeds 1 ms'
    )
