import numpy as np
import endf
import matplotlib.pyplot as plt

# Public API
__all__ = [
    "group_collapse_sigma",
    "unpack_covariance_ni_LB0",
    "unpack_covariance_ni_LB1",
    "unpack_covariance_ni_LB5",
    "unpack_covariance_ni_LB6",
    "unpack_covariance_ni"
]

def group_collapse_sigma(Ek, E_xs, sigma_xs):
    """
    Collapse fine-grid XS into group averages on {Ek}.
    Simple arithmetic average (could be flux-weighted in real use).
    """
    NI = len(Ek)-1
    sigma_g = np.zeros(NI)
    for i in range(NI):
        mask = (E_xs >= Ek[i]) & (E_xs < Ek[i+1])
        if np.any(mask):
            sigma_g[i] = np.mean(sigma_xs[mask])
        else:
            # fallback: interpolate midpoint
            midE = 0.5*(Ek[i] + Ek[i+1])
            sigma_g[i] = np.interp(midE, E_xs, sigma_xs)
    return sigma_g


def unpack_covariance_ni_LB0_old(mat, sub_section, ni, MT):
    data = mat.section_data[33, MT]
    block = data["subsections"][sub_section]["ni_subsections"][ni]

    Ek_cov = np.array(block["Ek"])
    NE = block.get("NE", len(Ek_cov))
    # NI = NE - 1
    dim = NE
    if "Fk" in block:
        Fdata = np.array(block["Fk"])
    else:
        raise NotImplementedError("No covariance numbers (Fkk) in block")

    C_abs = Fdata.reshape((dim, dim))

    return Ek_cov, Ek_cov, C_abs

def unpack_covariance_ni_LB0(mat, sub_section, ni, MT):
    data = mat.section_data[33, MT]
    block = data["subsections"][sub_section]["ni_subsections"][ni]

    Ek = np.array(block["Ek"])       # bin boundaries
    NE = block.get("NE", len(Ek))    # number of energy points
    Fk = np.array(block["Fk"])       # modal strengths (per bin)

    n_bins = len(Ek) - 1
    n_points = NE
    energies = Ek  # or construct a fine grid if needed

    # initialise covariance matrix
    C = np.zeros((n_points, n_points))

    for k in range(n_bins):
        # select which points fall in this energy bin
        mask_i = (energies >= Ek[k]) & (energies <= Ek[k+1])
        idx_i = np.where(mask_i)[0]

        # fill covariance for all pairs in the same bin
        for i in idx_i:
            for j in idx_i:
                C[i, j] += Fk[k]

    return Ek, Ek, C



def unpack_covariance_ni_LB1(mat, sub_section, ni, MT):
    data = mat.section_data[33, MT]
    block = data["subsections"][sub_section]["ni_subsections"][ni]
    MT1 = data["subsections"][sub_section]["MT1"]

    Ek_cov = np.array(block["Ek"])
    NE = block.get("NE", len(Ek_cov))
    NI = NE - 1

    # dim = NI
    # Get cross sections for MT and MT1
    xs_section = mat.section_data[3, MT]['sigma']
    E_xs, sigma_xs = xs_section.x, xs_section.y
    sigma_on_grid = group_collapse_sigma(Ek_cov, E_xs, sigma_xs)

    xs_section1 = mat.section_data[3, MT1]['sigma']
    E_xs1, sigma_xs1 = xs_section1.x, xs_section1.y
    sigma_on_grid1 = group_collapse_sigma(Ek_cov, E_xs1, sigma_xs1)

    # Get Fk data
    if "Fk" in block:
        Fdata = np.array(block["Fk"])
    else:
        raise NotImplementedError("No covariance numbers (Fk) in block")

    sigma = np.asarray(sigma_on_grid)
    sigma1 = np.asarray(sigma_on_grid1)

    if len(Ek_cov) != NI + 1:
        raise ValueError(f"Expected Ek of length NI+1 = {NI+1}, got {len(Ek_cov)}")
    if sigma.shape != (NE - 1,):
        raise ValueError(f"sigma must have shape (NE-1,), got {sigma.shape}")

    # Midpoints of collapsed groups
    group_edges = Ek_cov
    group_mids = 0.5 * (group_edges[:-1] + group_edges[1:])

    # Assign groups to covariance intervals
    interval_index = np.digitize(group_mids, Ek_cov) - 1  # values in 0..NI-1

    C_abs = np.zeros((NE-1, NE-1))
    for k, Fk in enumerate(Fdata):
        idx = np.where(interval_index == k)[0]
        for i in idx:
            for j in idx:
                C_abs[i, j] += Fk * sigma[i] * sigma1[j]

    return Ek_cov, Ek_cov, C_abs


def unpack_covariance_ni_LB5(mat, sub_section, ni, MT):
    data = mat.section_data[33, MT]

    block = data["subsections"][sub_section]["ni_subsections"][ni]
    # print(block)
    MT1 = data["subsections"][sub_section]["MT1"]

    Ek_cov = np.array(block["Ek"])
    NE = block.get("NE", len(Ek_cov))
    NI = NE - 1
    dim = NI

    if "LS" in block:
        LS = block["LS"]
    else:
        LS = 1

    xs_section = mat.section_data[3, MT]['sigma']
    E_xs, sigma_xs = xs_section.x, xs_section.y
    sigma_on_grid = group_collapse_sigma(Ek_cov, E_xs, sigma_xs)

    xs_section1 = mat.section_data[3, MT1]['sigma']
    E_xs1, sigma_xs1 = xs_section1.x, xs_section1.y
    sigma_on_grid1 = group_collapse_sigma(Ek_cov, E_xs1, sigma_xs1)

    # --- extract the covariance numbers ---
    if "Fkk" in block:
        Fdata = np.array(block["Fkk"])
    else:
        raise NotImplementedError("No covariance numbers (Fkk) in block")

    if LS == 0:  # full/asymmetric
        C = Fdata.reshape((dim, dim))
    elif LS == 1:  # symmetric packed
        C = np.zeros((dim, dim))
        idx = 0
        for i in range(dim):
            for j in range(i, dim):
                C[i, j] = Fdata[idx]
                C[j, i] = Fdata[idx]
                idx += 1
    else:
        raise ValueError("Unknown LS flag (must be 0 or 1)")

    # --- convert representation to absolute covariance ---
    sigma = np.asarray(sigma_on_grid)
    sigma1 = np.asarray(sigma_on_grid1)

    C_abs = np.outer(sigma, sigma1) * C

    return Ek_cov, Ek_cov, C_abs


def unpack_covariance_ni_LB6(mat, sub_section, ni, MT):
    """
    Unpack an LB=6 covariance block into an absolute covariance matrix.
    """
    data = mat.section_data[33, MT]
    block = data["subsections"][sub_section]["ni_subsections"][ni]
    MT1 = data["subsections"][sub_section]["MT1"]

    NER = block["NER"]
    NEC = block["NEC"]

    # Energies for rows (reaction MT) and columns (reaction MT1)
    Erow = np.array(block["ER"])
    Ecol = np.array(block["EC"])

    # Fractional covariance values, flattened row-major
    Fdata = np.array(block["Fkl"]).reshape((NER - 1, NEC - 1))

    # Cross sections on the same grids
    xs_section = mat.section_data[3, MT]['sigma']
    E_xs, sigma_xs = xs_section.x, xs_section.y
    sigma_row = group_collapse_sigma(Erow, E_xs, sigma_xs)

    xs_section1 = mat.section_data[3, MT1]['sigma']
    E_xs1, sigma_xs1 = xs_section1.x, xs_section1.y
    sigma_col = group_collapse_sigma(Ecol, E_xs1, sigma_xs1)

    sigma_row = np.asarray(sigma_row)
    sigma_col = np.asarray(sigma_col)

    # Absolute covariance: C_ij = Fkl_ij * σ_i * σ_j
    C_abs = np.zeros((NER - 1, NEC - 1))
    for i in range(NER - 1):
        for j in range(NEC - 1):
            C_abs[i, j] = Fdata[i, j] * sigma_row[i] * sigma_col[j]
    return Erow, Ecol, C_abs


def unpack_covariance_ni(mat, MT):
    results = []
    data = mat.section_data[33, MT]
    NL = data["NL"]

    for sub_section in range(0, NL):

        MT1 = data["subsections"][sub_section]["MT1"]
        if MT1 in range(850, 875) or (3, MT1) not in mat.sections:
            continue

        NI = data["subsections"][sub_section]["NI"]
        for ni in range(0, NI):
            print("gathering cov for MT", MT, " vs MT1", MT1 ," NI step = ", ni)
            # print(block)
            LB = data["subsections"][sub_section]["ni_subsections"][ni]["LB"]
            MT1 = data["subsections"][sub_section]["MT1"]

            # ----------------------------
            # LB = 0  (absolute covariance)
            # ----------------------------
            if LB == 0:
                Erow, Ecol, C_abs = unpack_covariance_ni_LB0(mat, sub_section, ni, MT)
            # ----------------------------
            # LB = 1  (relative covariance)
            # ----------------------------
            elif LB == 1:
                Erow, Ecol, C_abs = unpack_covariance_ni_LB1(mat, sub_section, ni, MT)
            # ----------------------------
            # LB = 5  (correlation matrix)
            # ----------------------------
            elif LB == 5:
                Erow, Ecol, C_abs = unpack_covariance_ni_LB5(mat, sub_section, ni, MT)
                # ----------------------------
                # LB = 6  (correlation matrix)
                # ----------------------------
            elif LB == 6:
                Erow, Ecol, C_abs = unpack_covariance_ni_LB6(mat, sub_section, ni, MT)
            # else:
            #     raise ValueError("Unknown LB ", LB)

            results.append({
                "Erow": Erow,
                "Ecol": Ecol,
                "C_abs": C_abs,
                "MT": MT,
                "MT1": MT1,
                "LB": LB,
                "NI_index": ni,
            })


    return results