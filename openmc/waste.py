from __future__ import annotations

import openmc
from openmc.data import half_life


def _waste_classification_nrc(mat: openmc.Material, metal: bool = True) -> str:
    """Classify a material according to 10 CFR 61.55 for near-surface disposal.

    Parameters
    ----------
    mat : openmc.Material
        The material to classify.
    metal : bool, optional
        Whether or not the material is in metal form. This changes the
        acceptable limits in Tables 1 and 2 for certain nuclides.

    Returns
    -------
    str
        The waste disposal classification, which can be "Class A", "Class B",
        "Class B", or "GTCC" (greater than class C).

    """
    # Specific activity limits from Table 1 in 10 CFR 61.55. For activated
    # metals, we manually add/modify the data.
    table1_limits = {
        'C14': (8.0, 'Ci/m3'),
        'Tc99': (3.0, 'Ci/m3'),
        'I129': (0.08, 'Ci/m3'),
        'Pu241': (3500.0, 'nCi/g'),
        'Cm242': (20000.0, 'nCi/g'),
        'Np237': (100.0, 'nCi/g'),
        'Pu238': (100.0, 'nCi/g'),
        'Pu239': (100.0, 'nCi/g'),
        'Pu240': (100.0, 'nCi/g'),
        'Pu242': (100.0, 'nCi/g'),
        'Pu244': (100.0, 'nCi/g'),
        'Am241': (100.0, 'nCi/g'),
        'Am243': (100.0, 'nCi/g'),
        'Cm243': (100.0, 'nCi/g'),
        'Cm244': (100.0, 'nCi/g'),
        'Cm245': (100.0, 'nCi/g'),
        'Cm246': (100.0, 'nCi/g'),
        'Cm247': (100.0, 'nCi/g'),
        'Cm248': (100.0, 'nCi/g'),
        'Bk247': (100.0, 'nCi/g'),
        'Cf249': (100.0, 'nCi/g'),
        'Cf250': (100.0, 'nCi/g'),
        'Cf251': (100.0, 'nCi/g'),
    }
    if metal:
        table1_limits['C14'] = (80.0, 'Ci/m3')
        table1_limits['Ni59'] = (220.0, 'Ci/m3')
        table1_limits['Nb94'] = (0.2, 'Ci/m3')

    # Specific activity limits from Table 2 in 10 CFR 61.55
    table2_limits = {
        'H3': (40.0, None, None),
        'Co60': (700.0, None, None),
        'Ni63': (3.5, 70.0, 700.0),
        'Sr90': (0.04, 150.0, 7000.0),
        'Cs137': (1.0, 44.0, 4600.0),
    }
    if metal:
        table2_limits['Ni63'] = (35.0, 700.0, 7000.0)

    # Add radionuclides with half-lives < 5 years to list for Table 2
    five_years = 60.0 * 60.0 * 24.0 * 365.25 * 5.0
    for nuc in mat.get_nuclides():
        if half_life(nuc) is not None and half_life(nuc) < five_years:
            table2_limits[nuc] = (700.0, None, None)

    # Determine which nuclides are present in Table 1 and Table 2
    table1_nuclides_present = any(x in table1_limits for x in mat.get_nuclides())
    table2_nuclides_present = any(x in table2_limits for x in mat.get_nuclides())

    # Determine metrics based on Tables 1 and 2 using sum of fractions rule for
    # mixture of radionuclides from §61.55(a)(7)
    g_cm3 = mat.get_mass_density()
    ratio1 = 0.0
    ratio2 = [0.0, 0.0, 0.0]
    for nuc, value in mat.get_activity(units="Ci/m3", by_nuclide=True).items():
        if nuc in table1_limits:
            limit, units = table1_limits[nuc]
            if units == 'nCi/g':
                value *= 1e9 / (1e6*g_cm3)  # Convert
            ratio1 += value / limit

        if nuc in table2_limits:
            for i, limit in enumerate(table2_limits[nuc]):
                if limit is not None:
                    ratio2[i] += value / limit

    # Helper function for classifying based on Table 2
    def classify_table2(col1, col2, col3):
        if col1 < 1.0:
            return "Class A"
        elif col2 < 1.0:
            return "Class B"
        elif col3 < 1.0:
            return "Class C"
        else:
            return "GTCC"

    if table1_nuclides_present and table2_nuclides_present:
        # Classification based on §61.55(a)(5)
        if ratio1 < 0.1:
            return classify_table2(*ratio2)
        elif ratio1 < 1.0:
            return "Class C" if ratio2[2] < 1.0 else "GTCC"
        else:
            return "GTCC"

    elif table1_nuclides_present:
        # Classification based on §61.55(a)(3)
        if ratio1 < 0.1:
            return "Class A"
        elif ratio1 < 1.0:
            return "Class C"
        else:
            return "GTCC"

    elif table2_nuclides_present:
        # Classification based on §61.55(a)(4)
        return classify_table2(*ratio2)

    else:
        # Classification based on §61.55(a)(6)
        return "Class A"


def _waste_classification_fetter(mat: openmc.Material) -> str:
    """Classify a material according to specific activity limits from Fetter.

    This function classifies a material based on the specific activity limits
    provided in the paper S. Fetter, E. T. Chang, and F. M. Mann, "Long-term
    radioactive waste from fusion reactors: Part II", Fusion Eng. Des., 13,
    239-246 (1990). https://doi.org/10.1016/0920-3796(90)90104-E. Because the
    NRC regulations in 10 CFR 61.55 do not cover all long-lived radionuclides
    relevant to fusion applications, this paper used the NRC methodology to
    calculate specific activity limits for an expanded set of radionuclides.

    Parameters
    ----------
    mat : openmc.Material
        The material to classify.

    Returns
    -------
    str
        The waste disposal classification, which can be either "Class C" or
        "GTCC" (greater than class C).

    """
    # Specific activity limits for class C disposal for all radionuclides with
    # half-lives between 5 years and 1e12 years from Table 2 in Fetter
    fetter_limit = {
        "Be10": 5.0e3,
        "C14": 6.0e2,
        "Al26": 9.0e-2,
        "Si32": 6.0e2,
        "Cl36": 1.0e1,
        "Ar39": 2.0e4,
        "Ar42": 2.0e4,
        "K40": 2.0e0,
        "Ca41": 1.0e4,
        "Ti44": 2.0e2,
        "Fe60": 1.0e-1,
        "Co60": 3.0e8,
        "Ni59": 9.0e2,
        "Ni63": 7.0e5,
        "Se79": 5.0e1,
        "Kr81": 3.0e1,
        "Sr90": 8.0e5,
        "Nb91": 2.0e2,
        "Nb92": 2.0e-1,
        "Nb94": 2.0e-1,
        "Mo93": 4.0e3,
        "Tc97": 4.0e-1,
        "Tc98": 1.0e-2,
        "Tc99": 6.0e-2,
        "Pb107": 9.0e2,
        "Ag108_m1": 3.0e0,
        "Sn121_m1": 7.0e5,
        "Sn126": 1.0e-1,
        "I129": 2.0e0,
        "Cs137": 5.0e4,
        "Ba133": 2.0e8,
        "La137": 2.0e2,
        "Sm151": 5.0e7,
        "Eu150_m1": 3.0e3,
        "Eu152": 3.0e5,
        "Eu154": 5.0e6,
        "Gd148": 2.0e5,
        "Gd150": 2.0e3,
        "Tb157": 5.0e3,
        "Tb158": 4.0e0,
        "Dy154": 1.0e3,
        "Ho166_m1": 2.0e-1,
        "Hf178_m1": 9.0e3,
        "Hf182": 2.0e-1,
        "Re186_m1": 2.0e1,
        "Ir192_m1": 1.0e0,
        "Pt193": 2.0e8,
        "Hg194": 5.0e-1,
        "Pb202": 6.0e-1,
        "Pb210": 3.0e7,
        "Bi207": 9.0e3,
        "Bi208": 8.0e-2,
        "Bi210_m1": 1.0e0,
        "Po209": 3.0e3,
        "Ra226": 1.0e-1,
        "Ra228": 3.0e7,
        "Ac227": 5.0e5,
        "Th229": 2.0e0,
        "Th230": 3.0e-1,
        "Th232": 1.0e-1,
        "Pa231": 7.0e-1,
        "U232": 3.0e1,
        "U233": 2.0e1,
        "U234": 9.0e1,
        "U235": 2.0e0,
        "Np236": 1.0e0,
        "Np237": 1.0e0,
        "Pu238": 7.0e4,
        "Pu239": 1.0e3,
        "Pu240": 1.0e3,
        "Pu241": 2.0e3,
        "Pu242": 1.0e3,
        "Pu244": 9.0e-1,
        "Am241": 5.0e1,
        "Am242_m1": 3.0e2,
        "Am243": 2.0e0,
        "Cm243": 6.0e2,
        "Cm244": 5.0e5,
        "Cm245": 5.0e0,
        "Cm246": 8.0e2,
        "Cm248": 8.0e2,
    }

    # Calculate the sum of the fractions of the activity of each radionuclide
    # compared to the Fetter limits
    ratio = 0.0
    for nuc, ci_m3 in mat.get_activity(units="Ci/m3", by_nuclide=True).items():
        if nuc in fetter_limit:
            ratio += ci_m3 / fetter_limit[nuc]

    # Classify based on the activity fraction
    return "Class C" if ratio < 1.0 else "GTCC"
