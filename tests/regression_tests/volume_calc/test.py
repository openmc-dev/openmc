import os
import glob

import numpy as np
import pytest

import openmc

from tests.testing_harness import PyAPITestHarness


class VolumeTest(PyAPITestHarness):

    def __init__(self, is_ce, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.exp_std_dev = 1e-01
        self.exp_rel_err = 1e-01
        self.exp_variance = 5e-02
        self.is_ce = is_ce
        if not is_ce:
            self.inputs_true = "inputs_true_mg.dat"

        # Define materials
        water = openmc.Material(1)
        water.add_nuclide("H1", 2.0)
        water.add_nuclide("O16", 1.0)
        water.add_nuclide("B10", 0.0001)
        if self.is_ce:
            water.add_s_alpha_beta("c_H_in_H2O")
        water.set_density("g/cc", 1.0)

        fuel = openmc.Material(2)
        fuel.add_nuclide("U235", 1.0)
        fuel.add_nuclide("Mo99", 0.1)
        fuel.set_density("g/cc", 4.5)

        materials = openmc.Materials((water, fuel))
        if not self.is_ce:
            materials.cross_sections = "mg_lib.h5"
        self._model.materials = materials

        cyl = openmc.ZCylinder(surface_id=1, r=1.0, boundary_type="vacuum")
        top_sphere = openmc.Sphere(surface_id=2, z0=5.0, r=1.0, boundary_type="vacuum")
        top_plane = openmc.ZPlane(surface_id=3, z0=5.0)
        bottom_sphere = openmc.Sphere(
            surface_id=4, z0=-5.0, r=1.0, boundary_type="vacuum"
        )
        bottom_plane = openmc.ZPlane(surface_id=5, z0=-5.0)

        # Define geometry
        inside_cyl = openmc.Cell(1, fill=fuel, region=-cyl & -top_plane & +bottom_plane)
        top_hemisphere = openmc.Cell(2, fill=water, region=-top_sphere & +top_plane)
        bottom_hemisphere = openmc.Cell(
            3, fill=water, region=-bottom_sphere & -top_plane
        )
        root = openmc.Universe(0, cells=(inside_cyl, top_hemisphere, bottom_hemisphere))

        self._model.geometry = openmc.Geometry(root)

        # Set up stochastic volume calculation
        ll, ur = root.bounding_box
        vol_calcs = [
            openmc.VolumeCalculation(list(root.cells.values()), 100000),
            openmc.VolumeCalculation([water, fuel], 100000, ll, ur),
            openmc.VolumeCalculation([root], 100000, ll, ur),
            openmc.VolumeCalculation(list(root.cells.values()), 100),
            openmc.VolumeCalculation([water, fuel], 100, ll, ur),
            openmc.VolumeCalculation(list(root.cells.values()), 100),
        ]

        vol_calcs[3].set_trigger(self.exp_std_dev, "std_dev")

        vol_calcs[4].set_trigger(self.exp_rel_err, "rel_err")

        vol_calcs[5].set_trigger(self.exp_variance, "variance")

        # Define settings
        settings = openmc.Settings()
        settings.run_mode = "volume"
        if not self.is_ce:
            settings.energy_mode = "multi-group"
        settings.volume_calculations = vol_calcs
        self._model.settings = settings

        # Create the MGXS file if necessary
        if not self.is_ce:
            groups = openmc.mgxs.EnergyGroups(group_edges=[0.0, 20.0e6])
            mg_xs_file = openmc.MGXSLibrary(groups)

            nu = [2.0]
            fiss = [1.0]
            capture = [1.0]
            absorption_fissile = np.add(fiss, capture)
            absorption_other = capture
            scatter = np.array([[[1.0]]])
            total_fissile = np.add(absorption_fissile, np.sum(scatter[:, :, 0], axis=1))
            total_other = np.add(absorption_other, np.sum(scatter[:, :, 0], axis=1))
            chi = [1.0]

            for iso in ["H1", "O16", "B10", "Mo99", "U235"]:
                mat = openmc.XSdata(iso, groups)
                mat.order = 0
                mat.atomic_weight_ratio = (
                    openmc.data.atomic_mass(iso) / openmc.data.NEUTRON_MASS
                )
                mat.set_scatter_matrix(scatter)
                if iso == "U235":
                    mat.set_nu_fission(np.multiply(nu, fiss))
                    mat.set_absorption(absorption_fissile)
                    mat.set_total(total_fissile)
                    mat.set_chi(chi)
                else:
                    mat.set_absorption(absorption_other)
                    mat.set_total(total_other)
                mg_xs_file.add_xsdata(mat)
            mg_xs_file.export_to_hdf5("mg_lib.h5")

    def _cleanup(self):
        super()._cleanup()
        output = ["mg_lib.h5"]
        for f in output:
            if os.path.exists(f):
                os.remove(f)

    def _get_results(self):
        outstr = ""
        for i, filename in enumerate(sorted(glob.glob("volume_*.h5"))):
            outstr += "Volume calculation {}\n".format(i)

            # Read volume calculation results
            volume_calc = openmc.VolumeCalculation.from_hdf5(filename)

            outstr += "Trigger Type: {}\n".format(volume_calc.trigger_type)
            outstr += "Trigger threshold: {}\n".format(volume_calc.threshold)
            outstr += "Iterations: {}\n".format(volume_calc.iterations)

            if i == 3:
                assert volume_calc.trigger_type == "std_dev"
                assert volume_calc.threshold == self.exp_std_dev
            elif i == 4:
                assert volume_calc.trigger_type == "rel_err"
                assert volume_calc.threshold == self.exp_rel_err
            elif i == 5:
                assert volume_calc.trigger_type == "variance"
                assert volume_calc.threshold == self.exp_variance
            else:
                assert volume_calc.trigger_type is None
                assert volume_calc.threshold is None
                assert volume_calc.iterations == 1

            # if a trigger is applied, make sure the calculation satisfies the trigger
            for vol in volume_calc.volumes.values():
                if volume_calc.trigger_type == "std_dev":
                    assert vol.std_dev <= self.exp_std_dev
                if volume_calc.trigger_type == "rel_err":
                    assert vol.std_dev / vol.nominal_value <= self.exp_rel_err
                if volume_calc.trigger_type == "variance":
                    assert vol.std_dev * vol.std_dev <= self.exp_variance

            # Write cell volumes and total # of atoms for each nuclide
            for uid, volume in sorted(volume_calc.volumes.items()):
                outstr += "Domain {}: {} cm^3\n".format(uid, volume)
            outstr += str(volume_calc.atoms_dataframe) + "\n"

        return outstr

    def _test_output_created(self):
        pass


@pytest.mark.parametrize("is_ce", [True, False])
def test_volume_calc(is_ce):
    harness = VolumeTest(is_ce, "", model=openmc.Model())
    harness.main()
