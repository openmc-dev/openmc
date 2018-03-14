.. _io_data_wmp:

=================================
Windowed Multipole Library Format
=================================

**/version** (*char[]*)
  The format version of the file.  The current version is "v0.2"

**/nuclide/**
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

        The residues are in the order: total, competitive if present,
        absorption, fission. Complex numbers are stored by forming a type with
        ":math:`r`" and ":math:`i`" identifiers, similar to how `h5py`_ does it.
    - **end_E** (*double*)
        Highest energy the windowed multipole part of the library is valid for.
    - **formalism** (*int*)
        The formalism of the underlying data. Uses the `ENDF-6`_ format
        formalism numbers.
        
        .. table:: Table of supported formalisms.
        
            +-------------+------------------+
            | Formalism   | Formalism number |
            +=============+==================+
            | MLBW        | 2                |
            +-------------+------------------+
            | Reich-Moore | 3                |
            +-------------+------------------+
        
    - **l_value** (*int[]*)
        The index for a corresponding pole. Equivalent to the :math:`l` quantum
        number of the resonance the pole comes from :math:`+1`.
    - **pseudo_K0RS** (*double[]*)
        :math:`l` dependent value of

        .. math::
            \sqrt{\frac{2 m_n}{\hbar}}\frac{AWR}{AWR + 1} r_{s,l}

        Where :math:`m_n` is mass of neutron, :math:`AWR` is the atomic weight
        ratio of the target to the neutron, and :math:`r_{s,l}` is the
        scattering radius for a given :math:`l`.
    - **spacing** (*double*)
        .. math::
            \frac{\sqrt{E_{max}}- \sqrt{E_{min}}}{n_w}

        Where :math:`E_{max}` is the maximum energy the windows go up to.  This
        is not equivalent to the maximum energy for which the windowed multipole
        data is valid for.  It is slightly higher to ensure an integer number of
        windows. :math:`E_{min}` is the minimum energy and equivalent to
        ``start_E``, and :math:`n_w` is the number of windows, given by
        ``windows``.
    - **sqrtAWR** (*double*)
        Square root of the atomic weight ratio.
    - **start_E** (*double*)
        Lowest energy the windowed multipole part of the library is valid for.
    - **w_start** (*int[]*)
        The pole to start from for each window.
    - **w_end** (*int[]*)
        The pole to end at for each window.

.. _h5py: http://docs.h5py.org/en/latest/
.. _ENDF-6: https://www.oecd-nea.org/dbdata/data/manual-endf/endf102.pdf
