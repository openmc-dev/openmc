.. _io_data_wmp:

=================================
Windowed Multipole Library Format
=================================

**/**

:Attributes: - **filetype** (*char[]*) -- String indicating the type of file
             - **version** (*int[2]*) -- Major and minor version of the data

**/<nuclide name>/**

:Datasets:

    - **broaden_poly** (*int[]*)
        If 1, Doppler broaden curve fit for window with corresponding index.
        If 0, do not.
    - **curvefit** (*double[][][]*)
        Curve fit coefficients. Indexed by (reaction type, coefficient index,
        window index).
    - **data** (*complex[][]*)
        Complex poles and residues. Each pole has a corresponding set of
        residues. For example, the :math:`i`-th pole and corresponding residues
        are stored as
        
        .. math::
            \text{data}[:,i] = [\text{pole},~\text{residue}_1,~\text{residue}_2,
            ~\ldots]

        The residues are in the order: scattering, absorption, fission. Complex
        numbers are stored by forming a type with ":math:`r`" and ":math:`i`"
        identifiers, similar to how `h5py`_ does it.
    - **E_max** (*double*)
        Highest energy the windowed multipole part of the library is valid for.
    - **E_min** (*double*)
        Lowest energy the windowed multipole part of the library is valid for.
    - **spacing** (*double*)
        .. math::
            \frac{\sqrt{E_{max}} - \sqrt{E_{min}}}{n_w}

        Where :math:`E_{max}` is the maximum energy the windows go up to.
        :math:`E_{min}` is the minimum energy, and :math:`n_w` is the number of
        windows, given by ``windows``.
    - **sqrtAWR** (*double*)
        Square root of the atomic weight ratio.
    - **windows** (*int[][]*)
        The poles to start from and end at for each window. windows[i, 0] and
        windows[i, 1] are, respectively, the indexes (1-based) of the first and
        last pole in window i.

.. _h5py: http://docs.h5py.org/en/latest/
