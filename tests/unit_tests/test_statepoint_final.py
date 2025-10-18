from pathlib import Path
import pytest

import openmc
import os


def test_statepoint_final(run_in_tmpdir):
    # Create a minimal model
    mat = openmc.Material()
    mat.add_nuclide('U235', 1.0)
    mat.set_density('g/cm3', 4.5)
    sphere = openmc.Sphere(r=10.0, boundary_type='vacuum')
    cell = openmc.Cell(fill=mat, region=-sphere)
    model = openmc.Model()
    model.geometry = openmc.Geometry([cell])
    model.settings.batches = 6
    model.settings.inactive = 2
    model.settings.particles = 50

    # Specify when statepoints should be written and trigger final write
    model.settings.statepoint = {'batches': [2, 4, 6]}

    # Locate built openmc executable in the workspace build directory. If it
    # doesn't exist, skip the test.
    repo_root = Path(__file__).resolve().parents[2]
    built_openmc = repo_root / 'build' / 'bin' / 'openmc'
    if not built_openmc.exists():
        pytest.skip(f"openmc executable not found at {built_openmc}; skip integration test")

    # If cross section environment variable isn't set, skip this test.
    if 'OPENMC_CROSS_SECTIONS' not in os.environ:
        pytest.skip("OPENMC_CROSS_SECTIONS not set; skip integration test that requires data libraries")

    # Run model and ensure that statepoints are created and final file exists
    last_sp = model.run(openmc_exec=str(built_openmc))

    # Check numbered statepoints using glob
    sp_glob = sorted([p.name for p in Path('.').glob('statepoint.*.h5')])
    if not sp_glob:
        # Provide diagnostics to help debugging when files are missing
        files = sorted([p.name for p in Path('.').iterdir()])
        raise AssertionError(f"No statepoint.*.h5 files found in {Path('.').resolve()}. Files present: {files}")

    # Expect at least 3 numbered statepoints to have been written
    assert len(sp_glob) >= 3, f"Expected at least 3 statepoint files, found: {sp_glob}"

    # Check deterministic final statepoint exists and was returned
    final = Path('statepoint.final.h5')
    assert final.is_file(), f"Expected statepoint.final.h5 to exist. Found: {sp_glob}"
    assert last_sp is not None, "Model.run should return a Path to the last statepoint"
    assert Path(last_sp).name == final.name, "Model.run should prefer statepoint.final.h5 when present"

    # Verify that statepoint.final.h5 mirrors the newest numbered statepoint
    import glob, h5py

    numbered = sorted(glob.glob('statepoint.*.h5'))
    # Remove statepoint.final.h5 from the numbered list if present
    numbered = [p for p in numbered if not p.endswith('statepoint.final.h5')]
    assert numbered, "No numbered statepoint files found to compare with final"

    # Pick the newest numbered statepoint by modification time
    newest = max(numbered, key=lambda p: Path(p).stat().st_mtime)

    # Read current_batch from both files and compare
    with h5py.File(newest, 'r') as f_num, h5py.File(final, 'r') as f_final:
        # current_batch may be at root or under run_info
        def _get_cb(fh):
            if 'current_batch' in fh:
                return fh['current_batch'][()]
            if 'run_info' in fh and 'current_batch' in fh['run_info']:
                # return fh['run_info/current_batch'][()]
                return fh['run_info']['current_batch'][()]
            return None

        cb_num = _get_cb(f_num)
        cb_final = _get_cb(f_final)

    assert cb_num is not None and cb_final is not None, "current_batch missing from one of the statepoint files"
    assert cb_num == cb_final, (
        f"statepoint.final.h5 current_batch ({cb_final}) does not match newest numbered "
        f"statepoint {Path(newest).name} current_batch ({cb_num})"
    )
