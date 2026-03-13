from typing import Callable
from warnings import warn

import openmc.lib


class _KeffSearchControl:
    """Controller for keff search during depletion calculations.

    This class performs keff searches to maintain a target keff by adjusting a
    model parameter through a provided function.

    Parameters
    ----------
    operator : openmc.deplete.Operator
        Depletion operator instance
    function : Callable
        Function that modifies the model based on a parameter value
    x0 : float
        Initial lower bound for the keff search
    x1 : float
        Initial upper bound for the keff search
    bracket : list[float]
        Absolute bracketing interval lower and upper. If the keff search
        solution lies off these limits the closest limit will be set as new
        result.
    **search_kwargs : dict, optional
        Additional keyword arguments to pass to :meth:`openmc.Model.keff_search`

    """
    def __init__(self, operator, function: Callable, x0: float, x1: float, bracket: list[float], **search_kwargs):
        if len(bracket) != 2:
            raise ValueError(f"bracket must have exactly 2 elements, got {len(bracket)}")
        if bracket[0] >= bracket[1]:
            raise ValueError(f"bracket[0] must be < bracket[1], got {bracket}")
        self.x0 = x0
        self.x1 = x1
        self.operator = operator
        self.function = function
        self.search_kwargs = search_kwargs
        self.search_kwargs['x_min'] = bracket[0]
        self.search_kwargs['x_max'] = bracket[1]

    def run(self, x):
        """Perform keff search and update the atom density vector.

        Parameters
        ----------
        x : list of numpy.ndarray
            Current atom density vector (atoms per material)

        Returns
        -------
        root : float
            Parameter value that achieves target keff
        """
        root = self._search_for_keff()
        self._update_vec(x)
        return root

    def _search_for_keff(self) -> float:
        """Perform the keff search using the model's keff_search method.

        Returns
        -------
        float
            Parameter value that achieves target keff

        Raises
        ------
        ValueError
            If the keff search fails to converge
        """
        with openmc.lib.TemporarySession(self.operator.model):
            # Only pass the first 3 required args plus explicitly provided kwargs
            result = self.operator.model.keff_search(
                self.function, self.x0, self.x1, **self.search_kwargs
            )
        if not result.converged:
            raise ValueError(
                f"Search for keff failed to converge. "
                f"Termination reason: {result.flag}"
            )

        root = result.root

        # Check if root is outside the bracket bounds and give a warning
        if root < self.search_kwargs['x_min']:
            warn(f"keff search result ({root:.6f}) is below the lower bracket "
                 f"bound ({self.search_kwargs['x_min']:.6f}).", UserWarning)
        elif root > self.search_kwargs['x_max']:
            warn(f"keff search result ({root:.6f}) is above the upper bracket "
                 f"bound ({self.search_kwargs['x_max']:.6f}).", UserWarning)

        # Restore the number of initial batches
        openmc.lib.settings.set_batches(self.operator.model.settings.batches)

        return root

    def _update_vec(self, x):
        """Update the atom density vector from openmc.lib.materials and AtomNumber object.

        The depletion vector ``x`` is rank-local, matching the materials owned
        by ``self.operator.number`` on the current MPI rank. We therefore only
        update entries for locally owned materials using the compositions
        currently stored in ``openmc.lib.materials``.

        Parameters
        ----------
        x : list of numpy.ndarray
            Atom density vector to update (atoms per material)

        """
        number = self.operator.number

        for mat_idx, mat in enumerate(number.materials):
            lib_material = openmc.lib.materials[int(mat)]
            nuclides = lib_material.nuclides
            densities = 1e24 * lib_material.densities
            volume = number.get_mat_volume(mat)

            for nuc_idx, nuc in enumerate(number.burnable_nuclides):
                if nuc in nuclides:
                    lib_nuc_idx = nuclides.index(nuc)
                    atom_density = densities[lib_nuc_idx]
                else:
                    atom_density = number.get_atom_density(mat, nuc)
                x[mat_idx][nuc_idx] = atom_density * volume
