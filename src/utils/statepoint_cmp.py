#!/usr/bin/env python

import sys

from numpy.testing import assert_allclose, assert_equal

from statepoint import StatePoint

if len(sys.argv) > 2:
    path1 = sys.argv[1]
    path2 = sys.argv[2]
else:
    raise

# Create StatePoint objects
sp1 = StatePoint(path1)
sp2 = StatePoint(path2)

# Read tally results
sp1.read_results()
sp2.read_results()

# Compare header information
assert sp1.revision == sp2.revision
assert sp1.version == sp2.version
assert sp1.seed == sp2.seed
assert sp1.run_mode == sp2.run_mode
assert sp1.n_particles == sp2.n_particles
assert sp1.n_batches == sp2.n_batches

if sp1.run_mode == 2:
    # Compare criticality information
    assert sp1.n_inactive == sp2.n_inactive
    assert sp1.gen_per_batch == sp2.gen_per_batch
    assert sp1.current_batch == sp2.current_batch

    # Compare keff results
    assert_allclose(sp1.k_batch, sp2.k_batch)

    # Compare entropy results
    assert_allclose(sp1.entropy, sp2.entropy)

# Compare global tallies
assert_allclose(sp1.global_tallies, sp2.global_tallies)

# Compare meshes
assert len(sp1.meshes) == len(sp2.meshes)
for m1, m2 in zip(sp1.meshes, sp2.meshes):
    assert m1.type == m2.type
    assert_equal(m1.dimension, m2.dimension)
    assert_equal(m1.lower_left, m2.lower_left)
    assert_equal(m1.upper_right, m2.upper_right)
    assert_equal(m1.width, m2.width)

# Compare tallies
assert len(sp1.tallies) == len(sp2.tallies)
for t1, t2 in zip(sp1.tallies, sp2.tallies):
    # Compare size of tallies
    assert t1.total_score_bins == t2.total_score_bins
    assert t1.total_filter_bins == t2.total_filter_bins

    # Compare filters
    assert len(t1.filters) == len(t2.filters)
    for f1, f2 in zip(t1.filters.values(), t2.filters.values()):
        assert f1.type == f2.type
        assert f1.length == f2.length
        assert_equal(f1.bins, f2.bins)

    # Compare nuclide and score bins
    assert_equal(t1.nuclides, t2.nuclides)
    assert_equal(t1.scores, t2.scores)

    # Compare tally results
    assert_allclose(t1.results, t2.results)

# If criticality, compare source sites
if sp1.run_mode == 2:
    sp1.read_source()
    sp2.read_source()

    assert len(sp1.source) == len(sp2.source)
    for s1, s2 in zip(sp1.source, sp2.source):
        assert s1.weight == s2.weight
        assert_allclose(s1.xyz, s2.xyz)
        assert_allclose(s1.uvw, s2.uvw)
        assert s1.E == s2.E
