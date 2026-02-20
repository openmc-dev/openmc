import openmc.lib
from openmc.mpi import comm
from typing import Callable
from warnings import warn

class KeffSearchControl:
    """Controller for keff search during depletion calculations.
    
    This class performs keff searches to maintain a target keff by adjusting
    a model parameter through a provided function.
    
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
        Absolute bracketing interval lower and upper
        if keff search solution lies off these limits the closest
        limit will be set as new result.
    search_kwargs : dict, optional
        Additional keyword arguments to pass to `model.keff_search`
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

    def search_for_keff(self, x, step_index):
        """Perform keff search and update the atom density vector.
        
        Parameters
        ----------
        x : list of numpy.ndarray
            Current atom density vector (atoms per material)
        step_index : int
            Current depletion step index
            
        Returns
        -------
        x : list of numpy.ndarray
            Updated atom density vector
        root : float
            Parameter value that achieves target keff
        """
        root = self._search_for_keff()
        x = self._update_vec(x)
        return x, root
        
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
        with openmc.lib.TemporarySession(model=self.operator.model) as session:
            # Only pass the first 3 required args plus explicitly provided kwargs
            result = self.operator.model.keff_search(
                self.function, 
                self.x0, 
                self.x1, 
                **self.search_kwargs
            )
        if not result.converged:
            raise ValueError(
                f"Search for keff failed to converge. "
                f"Termination reason: {result.termination_reason}"
            )
        
        root = result.root
        
        # Check if root is outside the bracket bounds and give a warning
        if root < self.search_kwargs['x_min']:
            warn(
                f"keff search result ({root:.6f}) is below the lower bracket bound "
                f"({self.search_kwargs['x_min']:.6f}). Clamping to bracket lower bound.",
                UserWarning
            )

        elif root > self.search_kwargs['x_max']:
            warn(
                f"keff search result ({root:.6f}) is above the upper bracket bound "
                f"({self.search_kwargs['x_max']:.6f}). Clamping to bracket upper bound.",
                UserWarning
            )
        
        #restore the number of initial batches
        openmc.lib.settings.set_batches(self.operator.model.settings.batches)

        return root

    def _update_vec(self, x):
        """Update the atom density vector from the operator's number object.
        
        This method synchronizes the atom densities across all MPI ranks by
        broadcasting the number object from each rank and updating the x vector
        with the current atom densities and volumes.
        
        Parameters
        ----------
        x : list of numpy.ndarray
            Atom density vector to update (atoms per material)
            
        Returns
        -------
        list of numpy.ndarray
            Updated atom density vector synchronized across all ranks
        """
        for rank in range(comm.size):
            number_i = comm.bcast(self.operator.number, root=rank)
            
            for mat_idx, mat in enumerate(number_i.materials):
                for nuc in number_i.nuclides:
                    if nuc in number_i.burnable_nuclides:
                        nuc_idx = number_i.burnable_nuclides.index(nuc)
                        atom_density = number_i.get_atom_density(mat, nuc)
                        volume = number_i.get_mat_volume(mat)
                        x[mat_idx][nuc_idx] = atom_density * volume
        
            x = comm.bcast(x, root=rank)
        return x