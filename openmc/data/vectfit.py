import numpy as np
import skrf
from skrf.vectorFitting import VectorFitting

def evaluate(s, poles, residues, polys=None):
    """
    Pure-NumPy version of vectfit.evaluate reproduced with scikit-rf conventions.
    Parameters
    ----------
    s : (Ns,) real or complex                     # sample points (typically j*ω)
    poles : (N,) complex
    residues : (Nv, N) complex
    polys : (Nv, Nc) real-valued coefficients, Nc=0/1/2 ...
    Returns
    -------
    f : (Nv, Ns) real (matches the C++ behaviour)
    """
    s   = np.asarray(s, dtype=complex)
    piv = 1.0 / (s[None, :] - poles[:, None])     # (N, Ns)

    out = (residues @ piv).real                  # rational part
    if polys is not None and polys.size:
        for m, c_m in enumerate(polys.T):        # add Σ c_m s**m
            out += c_m[:, None] * (s**m).real
    return out

def vectfit(f, s, poles_init, weight=None, n_polys=0,
            skip_pole=False, skip_res=False):
    """
    Replaces vectfit.vectfit() using pure-Python scikit-rf.
    Only the subset required by OpenMC (n_polys 0/1, no skip flags) is shown.
    """
    Nv, Ns = f.shape
    nports = int(np.sqrt(Nv))           # scikit-rf expects Nv == nports**2
    if nports**2 != Nv:
        raise ValueError("VectorFitting needs Nv to be a perfect square.")

    # 1-build an N-port Network object from the response matrix
    S = np.zeros((Ns, nports, nports), dtype=complex)
    for i in range(nports):
        for j in range(nports):
            S[:, i, j] = f[i*nports + j]          # row-major mapping

    freq = skrf.Frequency.from_f(s/(2*np.pi), unit='hz')   # s=jω → f
    nw = skrf.Network(frequency=freq, s=S, z0=50)        # 50 Ω dummy

    # 2-run vector fitting
    vf = VectorFitting(nw)
    vf.poles = poles_init
    vf.vector_fit(init_pole_spacing='custom',
                  n_poles_real = np.count_nonzero(poles_init.imag == 0),
                  n_poles_cmplx= np.count_nonzero(poles_init.imag > 0))

    # 3-package results so OpenMC sees the same tuple
    fitted = np.vstack([vf.get_model_response(i//nports, i % nports)
                        for i in range(Nv)])
    rmserr = np.linalg.norm(fitted - f) / np.sqrt(f.size)

    # scikit-rf holds only degree-0/1 terms; stack them to look like polys
    polys  = np.stack([vf.constant_coeff,
                       vf.proportional_coeff], axis=1)[:, :n_polys]

    return vf.poles, vf.residues, polys, fitted, rmserr
