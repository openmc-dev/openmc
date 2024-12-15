import pytest

# Test configuration options for regression tests
config = {
    "event": False,
    "exe": "openmc",
    "mpi": False,
    "mpiexec": "mpiexec",
    "mpi_np": "2",
    "update": False,
    "build_inputs": False,
}


def assert_same_mats(res_ref, res_test):
    for mat in res_ref[0].index_mat:
        assert mat in res_test[0].index_mat, f"Material {mat} not in new results."
    for nuc in res_ref[0].index_nuc:
        assert nuc in res_test[0].index_nuc, f"Nuclide {nuc} not in new results."
    for mat in res_test[0].index_mat:
        assert mat in res_ref[0].index_mat, f"Material {mat} not in old results."
    for nuc in res_test[0].index_nuc:
        assert nuc in res_ref[0].index_nuc, f"Nuclide {nuc} not in old results."


def assert_atoms_equal(res_ref, res_test, tol=1e-5):
    for mat in res_test[0].index_mat:
        for nuc in res_test[0].index_nuc:
            _, y_test = res_test.get_atoms(mat, nuc)
            _, y_ref = res_ref.get_atoms(mat, nuc)
            assert y_test == pytest.approx(y_ref, rel=tol), (
                f"Atoms not equal for material {mat}, nuclide {nuc}\n"
                f"y_ref={y_ref}\ny_test={y_test}"
            )


def assert_reaction_rates_equal(res_ref, res_test, tol=1e-5):
    for reactions in res_test[0].rates:
        for mat in reactions.index_mat:
            for nuc in reactions.index_nuc:
                for rx in reactions.index_rx:
                    y_test = res_test.get_reaction_rate(mat, nuc, rx)[1]
                    y_ref = res_ref.get_reaction_rate(mat, nuc, rx)[1]
                    assert y_test == pytest.approx(y_ref, rel=tol), (
                        f"Reaction rate not equal for material {mat}, nuclide "
                        f"{nuc}, {rx}\ny_ref={y_ref}\ny_test={y_test}"
                    )
